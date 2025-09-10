"""
make_fel.py
-----------
Peptide-based Free Energy Landscape (FEL) utilities for oligomerization analysis.

Modes
=====
A) 1D FEL: F(n) from existing cluster CSVs (no trajectory needed)
   Input: cluster_results/cluster_sizes_per_frame.csv
   Output: fel_F_of_n.(csv|png|pdf)

B) 2D FEL: F(n, Rg) where Rg is the radius of gyration of the *largest* cluster per frame
   (requires topology+trajectory; recomputes clusters peptide-wise).
   Output: fel_F_of_n_Rg.(npy|png|pdf) and a CSV of bin centers + free energies.

Usage
=====
# Mode A (1D FEL from CSV; fast)
python make_clust_fel.py --indir results_cluster --temperature 300

# Mode B (2D FEL with Rg; recompute clusters)
python make_fel.py \
  --top topol.tpr --traj traj.xtc --groupby auto \
  --selection "protein and name CA" \
  --cutoff 0.6 --min_contacts 3 --stride 10 --pbc \
  --temperature 300 \
  --fel2d --outdir fel_results

Notes
=====
- FEL computed as F = -k_B T ln P + C (C chosen so min(F)=0).
- Add a small pseudocount to avoid log(0).
- Units: k_B in kJ/(mol*K), so F is in kJ/mol.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kB = 0.008314462618  # kJ/(mol*K)

def parse_args():
    p = argparse.ArgumentParser(description="Compute peptide-based FELs over oligomer size (and optional Rg).")
    p.add_argument("--indir", default="cluster_results", help="Directory with cluster CSVs (Mode A)")
    p.add_argument("--temperature", type=float, default=300.0, help="Temperature in K (default 300)")
    p.add_argument("--outdir", default=None, help="Output directory (defaults to --indir)")

    # Mode B (trajectory-based 2D FEL)
    p.add_argument("--fel2d", action="store_true", help="Enable 2D FEL over (n, Rg_largest) [requires traj]")
    p.add_argument("--top", help="Topology file (e.g., topol.tpr, PDB)")
    p.add_argument("--traj", help="Trajectory file (e.g., XTC, DCD)")
    p.add_argument("--selection", default="protein and name CA",
                   help='Atom selection for cluster detection and Rg (default: "protein and name CA")')
    p.add_argument("--groupby", choices=["auto","segid","chainID","molnum"], default="auto",
                   help="Grouping key for peptides (default: auto)")
    p.add_argument("--cutoff", type=float, default=0.6, help="Contact cutoff (units match trajectory coords)")
    p.add_argument("--min_contacts", type=int, default=3, help="Min atom-atom contacts to consider connected")
    p.add_argument("--stride", type=int, default=10, help="Use every Nth frame")
    p.add_argument("--pbc", action="store_true", help="Use PBC minimum-image distances")
    p.add_argument("--start_ns", type=float, default=None, help="Optional start time (ns)")
    p.add_argument("--stop_ns", type=float, default=None, help="Optional stop time (ns)")
    p.add_argument("--bins_n", type=int, default=9, help="Bins for n (default up to 9)")
    p.add_argument("--bins_Rg", type=int, default=40, help="Bins for Rg (default 40)")
    p.add_argument("--Rg_unit", choices=["nm","angstrom"], default="nm", help="Rg output unit (default nm)")
    p.add_argument("--pseudocount", type=float, default=1e-6, help="Probability pseudocount")
    return p.parse_args()

# ------------------ Mode A: FEL over n from CSV ------------------

def fel_1d_from_csv(indir, T, outdir=None, pseudocount=1e-6):
    indir = Path(indir)
    outdir = Path(outdir) if outdir else indir
    outdir.mkdir(parents=True, exist_ok=True)
    csv = indir / "cluster_sizes_per_frame.csv"
    if not csv.exists():
        raise FileNotFoundError(f"Missing {csv}. Run cluster_analysis.py first.")

    df = pd.read_csv(csv)
    if "largest_cluster" not in df.columns:
        raise ValueError("cluster_sizes_per_frame.csv is missing 'largest_cluster'. Update your analysis script.")

    n_vals = df["largest_cluster"].values.astype(int)
    # Build histogram over n >=1
    n_max = int(max(n_vals)) if len(n_vals)>0 else 1
    bins = np.arange(0.5, n_max+1.5, 1.0)  # centers at 1,2,...,n_max
    hist, edges = np.histogram(n_vals, bins=bins, density=False)
    counts = hist.astype(float)
    P = (counts + pseudocount) / (counts.sum() + pseudocount*len(counts))

    F = -kB * T * np.log(P)
    F -= np.nanmin(F)  # set minimum to 0
    centers = np.arange(1, n_max+1, 1)
    out = pd.DataFrame({"n": centers, "prob": P, "F_kJmol": F})
    out.to_csv(outdir / "fel_F_of_n.csv", index=False)

    # Plot
    plt.figure(figsize=(6,4))
    plt.plot(centers, F, marker="o")
    plt.xlabel("Largest cluster size, n")
    plt.ylabel("F(n) [kJ/mol]")
    plt.title(f"FEL over n (T={T:.0f} K)")
    plt.tight_layout()
    plt.savefig(outdir / "fel_F_of_n.png", dpi=300)
    plt.savefig(outdir / "fel_F_of_n.pdf")
    plt.close()
    print("Wrote 1D FEL to:", outdir)

# ------------------ Mode B helpers: clustering + Rg ---------------

def _unique_nonempty(arr):
    vals = []
    for v in arr:
        if isinstance(v, str):
            if v.strip() == "":
                continue
        vals.append(v)
    # preserve order
    seen = set()
    keep = []
    for v in vals:
        if v in seen:
            continue
        seen.add(v)
        keep.append(v)
    return keep

def _build_groups(u, selection, groupby):
    sel = u.select_atoms(selection)
    if sel.n_atoms == 0:
        raise ValueError(f"No atoms matched selection: {selection}")

    def groups_from(label, keys):
        gs = []
        used = set()
        for key in keys:
            if key in used:
                continue
            used.add(key)
            q = f"{label} '{key}'" if isinstance(key, str) else f"{label} {key}"
            ag = sel.select_atoms(q)
            if ag.n_atoms > 0:
                gs.append(ag)
        return gs

    # explicit choice
    if groupby in ("segid","chainID","molnum"):
        keys = getattr(sel.atoms, groupby+"s")
        keys = _unique_nonempty(keys)
        gs = groups_from(groupby, keys)
        if len(gs) >= 2:
            return gs
        raise ValueError(f"Found {len(gs)} group(s) with {groupby}")

    # auto: segid -> chainID -> molnum -> residue reset
    try:
        keys = _unique_nonempty(sel.atoms.segids)
        gs = groups_from("segid", keys)
        if len(gs) >= 2:
            return gs
    except Exception:
        pass
    try:
        keys = _unique_nonempty(sel.atoms.chainIDs)
        gs = groups_from("chainID", keys)
        if len(gs) >= 2:
            return gs
    except Exception:
        pass
    try:
        keys = _unique_nonempty(sel.atoms.molnums)
        gs = groups_from("molnum", keys)
        if len(gs) >= 2:
            return gs
    except Exception:
        pass

    # residue reset (fallback): split CA atoms where resid drops significantly
    ca = sel.select_atoms("name CA") if sel.select_atoms("name CA").n_atoms>0 else sel
    resids = ca.resids
    breaks = [0]
    for i in range(1, len(resids)):
        if resids[i] <= resids[i-1] - 10:
            breaks.append(i)
    breaks.append(len(ca))
    gs = []
    for b0, b1 in zip(breaks[:-1], breaks[1:]):
        idxs = ca.atoms.indices[b0:b1]
        ag = sel.atoms[np.isin(sel.atoms.indices, idxs)]
        if ag.n_atoms>0:
            gs.append(ag)
    if len(gs) >= 2:
        return gs
    raise ValueError("Could not determine peptide groups automatically.")

def _capped_pairs(ri, rj, cutoff, box):
    from MDAnalysis.lib.distances import capped_distance
    pairs = capped_distance(ri, rj, max_cutoff=cutoff, box=box, return_distances=False)
    return pairs.shape[0] if hasattr(pairs, "shape") else len(pairs)

def _connected_components(n, edges):
    parent = list(range(n))
    rank = [0]*n
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]]
            x=parent[x]
        return x
    def union(x,y):
        rx,ry=find(x),find(y)
        if rx==ry: return
        if rank[rx]<rank[ry]:
            parent[rx]=ry
        elif rank[rx]>rank[ry]:
            parent[ry]=rx
        else:
            parent[ry]=rx
            rank[rx]+=1
    for i,j in edges: union(i,j)
    groups = {}
    for i in range(n):
        r=find(i); groups.setdefault(r, []).append(i)
    return list(groups.values())

def _rg(ag, box=None):
    # Simple radius of gyration around the center of mass; ignore masses (uniform).
    coords = ag.positions
    com = coords.mean(axis=0)
    diffs = coords - com
    if box is not None:
        # unwrap minimum-image relative to COM for triclinic boxes
        # approximate: MDAnalysis has triclinic unwrap in distances; here we skip due to complexity.
        pass
    rg2 = (diffs*diffs).sum(axis=1).mean()
    return float(np.sqrt(rg2))

def fel_2d_recompute(top, traj, selection, groupby, cutoff, min_contacts, stride, pbc, T,
                     start_ns=None, stop_ns=None, bins_n=9, bins_Rg=40, Rg_unit="nm",
                     outdir=None, pseudocount=1e-6):
    try:
        import MDAnalysis as mda
    except Exception as e:
        print("MDAnalysis import failed. pip install MDAnalysis", file=sys.stderr)
        raise

    outdir = Path(outdir) if outdir else Path("fel_results")
    outdir.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(top, traj)
    groups = _build_groups(u, selection, groupby)
    n_groups = len(groups)

    # Determine cutoff units (heuristic like before): assume traj coords in nm if box length < 50
    cutoff_internal = cutoff
    inferred = None
    for ts in u.trajectory[:: max(1, len(u.trajectory)//5) ]:
        if ts.dimensions is not None and ts.dimensions[0] > 0:
            L = float(ts.dimensions[0])
            inferred = "angstrom" if L > 50 else "nm"
            break
    # If coords are Angstroms but cutoff provided in nm, scale
    if inferred == "angstrom":
        cutoff_internal = cutoff * 10.0

    Ns = []
    Rgs = []

    start_ps = start_ns*1000.0 if start_ns is not None else None
    stop_ps  = stop_ns*1000.0 if stop_ns  is not None else None

    for iframe, ts in enumerate(u.trajectory[::stride]):
        t_ps = float(ts.time) if hasattr(ts, "time") else None
        if start_ps is not None and t_ps is not None and t_ps < start_ps:
            continue
        if stop_ps is not None and t_ps is not None and t_ps > stop_ps:
            break

        # Build edges
        edges = []
        box = ts.dimensions if pbc else None
        for i in range(n_groups-1):
            gi = groups[i]
            for j in range(i+1, n_groups):
                gj = groups[j]
                c = _capped_pairs(gi.positions, gj.positions, cutoff_internal, box)
                if c >= min_contacts:
                    edges.append((i,j))

        comps = _connected_components(n_groups, edges)
        sizes = [len(c) for c in comps]
        if not sizes:
            Ns.append(0)
            Rgs.append(np.nan)
            continue
        # largest component
        imax = int(np.argmax(sizes))
        largest_nodes = comps[imax]
        ag = groups[largest_nodes[0]]
        for idx in largest_nodes[1:]:
            ag = ag + groups[idx]
        rg = _rg(ag, box=box)
        if inferred == "angstrom" and Rg_unit=="nm":
            rg = rg / 10.0
        Rgs.append(float(rg))
        Ns.append(int(sizes[imax]))

    # 2D histogram over (n, Rg)
    n_max = max(Ns) if Ns else 1
    n_bins = np.arange(0.5, n_max+1.5, 1.0)  # integer bins
    if isinstance(bins_n, int):
        # ensure we cover full integer range regardless of bins_n
        n_bins = np.arange(0.5, n_max+1.5, 1.0)
    rg_min = np.nanmin(Rgs) if len(Rgs)>0 else 0.0
    rg_max = np.nanmax(Rgs) if len(Rgs)>0 else 1.0
    # safe margins
    rg_pad = 0.05*(rg_max - rg_min + 1e-12)
    rg_edges = np.linspace(rg_min - rg_pad, rg_max + rg_pad, bins_Rg+1)

    # Mask NaNs
    Ns_arr = np.array(Ns)
    Rg_arr = np.array(Rgs)
    mask = np.isfinite(Ns_arr) & np.isfinite(Rg_arr) & (Ns_arr>0)
    H, n_edges, rg_edges = np.histogram2d(Ns_arr[mask], Rg_arr[mask], bins=[n_bins, rg_edges])
    H = H.astype(float)
    P = (H + pseudocount) / (H.sum() + pseudocount*H.size)
    F = -kB * T * np.log(P)
    F -= np.nanmin(F)

    # Save CSV with bin centers
    n_centers = 0.5*(n_edges[:-1] + n_edges[1:])
    rg_centers = 0.5*(rg_edges[:-1] + rg_edges[1:])
    # Flatten grid to CSV
    rows = []
    for i, n_c in enumerate(n_centers):
        for j, rg_c in enumerate(rg_centers):
            rows.append({"n_center": n_c, f"Rg_{Rg_unit}": rg_c, "F_kJmol": F[i,j]})
    fel2d_df = pd.DataFrame(rows)
    out_csv = outdir / "fel_F_of_n_Rg.csv"
    fel2d_df.to_csv(out_csv, index=False)

    # Save numpy array too
    np.save(outdir / "fel_F_of_n_Rg.npy", F)

    # Plot as heatmap
    plt.figure(figsize=(6,5))
    extent = [rg_edges[0], rg_edges[-1], n_edges[0], n_edges[-1]]
    # imshow expects [xmin,xmax,ymin,ymax]; origin lower to have small n at bottom
    plt.imshow(F.T, origin="lower", aspect="auto", extent=extent)
    plt.xlabel(f"Rg (largest cluster) [{Rg_unit}]")
    plt.ylabel("Cluster size, n (largest)")
    plt.title(f"F(n, Rg) at T={T:.0f} K")
    cbar = plt.colorbar()
    cbar.set_label("F [kJ/mol]")
    plt.tight_layout()
    plt.savefig(outdir / "fel_F_of_n_Rg.png", dpi=300)
    plt.savefig(outdir / "fel_F_of_n_Rg.pdf")
    plt.close()

    print("Wrote 2D FEL to:", outdir)

def main():
    args = parse_args()
    outdir = args.outdir or args.indir
    if args.fel2d:
        if not (args.top and args.traj):
            raise SystemExit("For --fel2d you must provide --top and --traj")
        fel_2d_recompute(args.top, args.traj, args.selection, args.groupby,
                         args.cutoff, args.min_contacts, args.stride, args.pbc,
                         args.temperature, args.start_ns, args.stop_ns,
                         args.bins_n, args.bins_Rg, args.Rg_unit,
                         outdir, args.pseudocount)
    else:
        fel_1d_from_csv(args.indir, args.temperature, outdir, args.pseudocount)

if __name__ == "__main__":
    main()
