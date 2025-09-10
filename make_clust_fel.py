#!/usr/bin/env python3
"""
make_fel.py  (v1.1)
-------------------
Peptide-based FEL utilities with optional monomer chemical-potential correction.

New:
- --mu_mode none|fraction|concentration
  * fraction: F_corr(n) = F_sim(n) + n*kT*ln((N-n)/N)
  * concentration: mu = kT*ln(c1/c0), c1=(N-n)/(N_A * V_L);  F_corr(n)=F_sim(n) - n*mu
    Requires --N_total and either --box_nm3 or --volume_from_traj.

- --N_total, --box_nm3, --std_conc_M (default 1.0), --volume_from_traj

Modes remain:
A) 1D FEL over n from cluster_sizes_per_frame.csv
B) 2D FEL over (n, Rg) with --fel2d (no mu correction applied to 2D map by default).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#kB = 0.008314462618  # kJ/(mol*K)
kB = 0.001985875  # kcal/(mol*K)
N_A = 6.02214076e23  # 1/mol
NM3_TO_L = 1e-24     # 1 nm^3 = 1e-24 L

def parse_args():
    p = argparse.ArgumentParser(description="Compute peptide-based FELs over oligomer size (and optional Rg).")
    p.add_argument("--indir", default="cluster_results", help="Directory with cluster CSVs (Mode A)")
    p.add_argument("--temperature", type=float, default=300.0, help="Temperature in K (default 300)")
    p.add_argument("--outdir", default=None, help="Output directory (defaults to --indir)")

    # Chemical potential correction options (1D FEL only)
    p.add_argument("--mu_mode", choices=["none","fraction","concentration"], default="none",
                   help="Apply monomer chemical-potential correction to 1D FEL (default: none)")
    p.add_argument("--N_total", type=int, default=None, help="Total number of peptides N (required for mu correction)")
    p.add_argument("--box_nm3", type=float, default=None, help="Simulation box volume in nm^3 (for mu=concentration)")
    p.add_argument("--std_conc_M", type=float, default=1.0, help="Standard-state concentration c0 in M (default 1.0)")
    p.add_argument("--volume_from_traj", action="store_true",
                   help="Estimate average box volume from trajectory (requires --top/--traj)")

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

def _estimate_volume_from_traj(top, traj, stride=50):
    try:
        import MDAnalysis as mda
    except Exception as e:
        raise RuntimeError("Need MDAnalysis installed for --volume_from_traj") from e
    u = mda.Universe(top, traj)
    vols = []
    for ts in u.trajectory[::stride]:
        # ts.dimensions[:3] are triclinic vectors' lengths; volume = a*b*c*sqrt(1 + 2cosαcosβcosγ - cos^2α - cos^2β - cos^2γ)
        # MDAnalysis provides volume as ts.volume (nm^3) when available
        if getattr(ts, "volume", None) is not None:
            vols.append(float(ts.volume))
        else:
            # Fallback to orthorhombic approximation
            a, b, c = ts.dimensions[:3]
            vols.append(float(a*b*c))
    if not vols:
        raise RuntimeError("Could not read volume from trajectory.")
    return float(np.mean(vols))

def fel_1d_from_csv(indir, T, outdir=None, pseudocount=1e-6,
                    mu_mode="none", N_total=None, box_nm3=None, std_conc_M=1.0,
                    volume_from_traj=False, top=None, traj=None):
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
    n_max = int(max(n_vals)) if len(n_vals)>0 else 1
    centers = np.arange(1, n_max+1, 1)
    bins = np.arange(0.5, n_max+1.5, 1.0)
    hist, _ = np.histogram(n_vals, bins=bins, density=False)
    counts = hist.astype(float)
    P = (counts + pseudocount) / (counts.sum() + pseudocount*len(counts))
    Fsim = -kB * T * np.log(P)
    Fsim -= np.nanmin(Fsim)

    # --- Chemical potential correction (optional) ---
    Fcorr = Fsim.copy()
    extra_note = ""
    if mu_mode != "none":
        if N_total is None:
            raise ValueError("--N_total is required for --mu_mode fraction|concentration")

        if mu_mode == "fraction":
            # F_corr(n) = F_sim(n) + n*kT*ln((N-n)/N)
            frac = (N_total - centers) / float(N_total)
            # avoid log(0) at n=N
            frac = np.clip(frac, 1e-12, None)
            Fcorr = Fsim + kB*T * (centers * np.log(frac))
            extra_note = "Applied fraction-based mu correction: F_corr(n) = F_sim(n) + n kT ln((N-n)/N)."

        elif mu_mode == "concentration":
            # Need a volume in liters
            if volume_from_traj:
                if top is None or traj is None:
                    raise ValueError("--volume_from_traj requires --top and --traj")
                V_nm3 = _estimate_volume_from_traj(top, traj)
            else:
                if box_nm3 is None:
                    raise ValueError("--box_nm3 (nm^3) or --volume_from_traj must be provided for mu=concentration")
                V_nm3 = box_nm3
            V_L = V_nm3 * NM3_TO_L
            # c1(n) in mol/L
            c1 = (N_total - centers) / (N_A * V_L)
            c1 = np.clip(c1, 1e-30, None)  # avoid log(0)
            mu = kB*T * np.log(c1 / float(std_conc_M))
            Fcorr = Fsim - centers * mu
            extra_note = f"Applied concentration-based mu correction with N={N_total}, V={V_nm3:.3e} nm^3, c0={std_conc_M} M."

        # Shift so min=0 for readability
        Fcorr -= np.nanmin(Fcorr)

    out = pd.DataFrame({"n": centers, "prob": P, "F_sim_kJmol": Fsim, "F_corr_kJmol": Fcorr})
    out.to_csv(outdir / "fel_F_of_n.csv", index=False)

    # Plot both (if corrected, show both curves)
    plt.figure(figsize=(6,4))
    lbl_sim = "F_sim(n)"
    plt.plot(centers, Fsim, marker="o", linestyle="-", label=lbl_sim)
    if mu_mode != "none":
        plt.plot(centers, Fcorr, marker="s", linestyle="--", label="F_corr(n)")
        plt.legend()
    plt.xlabel("Largest cluster size, n")
    plt.ylabel("F [kcal/mol]")
    title = f"FEL over n (T={T:.0f} K)"
    if mu_mode != "none":
        title += f"  [{mu_mode} μ-correction]"
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outdir / "fel_F_of_n.png", dpi=300)
    plt.savefig(outdir / "fel_F_of_n.pdf")
    plt.close()

    # Write a small README note
    with open(outdir / "fel_README.txt", "w") as f:
        f.write("FEL over n computed as F = -kT ln P + C; minimum set to 0.\n")
        if mu_mode != "none":
            f.write(extra_note + "\n")

# ------------------ Mode B: 2D FEL remains the same ---------------

# (Import the previous v1.0 2D code by reading from the same file, to avoid duplication here.)
# For brevity in this snippet, we just re-define the same functions from the prior version.

def _unique_nonempty(arr):
    vals = []
    for v in arr:
        if isinstance(v, str):
            if v.strip() == "":
                continue
        vals.append(v)
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

    if groupby in ("segid","chainID","molnum"):
        keys = getattr(sel.atoms, groupby+"s")
        keys = _unique_nonempty(keys)
        gs = groups_from(groupby, keys)
        if len(gs) >= 2:
            return gs
        raise ValueError(f"Found {len(gs)} group(s) with {groupby}")

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
    coords = ag.positions
    com = coords.mean(axis=0)
    diffs = coords - com
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

    cutoff_internal = cutoff
    inferred = None
    for ts in u.trajectory[:: max(1, len(u.trajectory)//5) ]:
        if ts.dimensions is not None and ts.dimensions[0] > 0:
            L = float(ts.dimensions[0])
            inferred = "angstrom" if L > 50 else "nm"
            break
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
            Ns.append(0); Rgs.append(np.nan); continue
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

    n_max = max(Ns) if Ns else 1
    n_bins = np.arange(0.5, n_max+1.5, 1.0)
    rg_min = np.nanmin(Rgs) if len(Rgs)>0 else 0.0
    rg_max = np.nanmax(Rgs) if len(Rgs)>0 else 1.0
    rg_pad = 0.05*(rg_max - rg_min + 1e-12)
    rg_edges = np.linspace(rg_min - rg_pad, rg_max + rg_pad, bins_Rg+1)

    Ns_arr = np.array(Ns)
    Rg_arr = np.array(Rgs)
    mask = np.isfinite(Ns_arr) & np.isfinite(Rg_arr) & (Ns_arr>0)
    H, n_edges, rg_edges = np.histogram2d(Ns_arr[mask], Rg_arr[mask], bins=[n_bins, rg_edges])
    H = H.astype(float)
    P = (H + pseudocount) / (H.sum() + pseudocount*H.size)
    F = -kB * T * np.log(P)
    F -= np.nanmin(F)

    n_centers = 0.5*(n_edges[:-1] + n_edges[1:])
    rg_centers = 0.5*(rg_edges[:-1] + rg_edges[1:])
    rows = []
    for i, n_c in enumerate(n_centers):
        for j, rg_c in enumerate(rg_centers):
            rows.append({"n_center": n_c, f"Rg_{Rg_unit}": rg_c, "F_kJmol": F[i,j]})
    fel2d_df = pd.DataFrame(rows)
    out_csv = outdir / "fel_F_of_n_Rg.csv"
    fel2d_df.to_csv(out_csv, index=False)

    np.save(outdir / "fel_F_of_n_Rg.npy", F)

    plt.figure(figsize=(6,5))
    extent = [rg_edges[0], rg_edges[-1], n_edges[0], n_edges[-1]]
    plt.imshow(F.T, origin="lower", aspect="auto", extent=extent)
    plt.xlabel(f"Rg (largest cluster) [{Rg_unit}]")
    plt.ylabel("Cluster size, n (largest)")
    plt.title(f"F(n, Rg) at T={T:.0f} K")
    cbar = plt.colorbar(); cbar.set_label("F [kcal/mol]")
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
        fel_1d_from_csv(args.indir, args.temperature, outdir, args.pseudocount,
                        mu_mode=args.mu_mode, N_total=args.N_total, box_nm3=args.box_nm3,
                        std_conc_M=args.std_conc_M, volume_from_traj=args.volume_from_traj,
                        top=args.top, traj=args.traj)

if __name__ == "__main__":
    main()

