import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import math
from itertools import combinations
import networkx as nx
from tqdm import tqdm
import argparse
import os

# target = 1
cutoff_ca = 8.5
# min_res_pairs = 3 if target <= 5 else 4

# --- choose per-n settings (keep your hexamer defaults) ---
def params_for(n):
    # alpha = fraction of complete-graph edges required
    if n <= 5:
        return dict(alpha=0.55, min_k=2, max_diam=None)  # looser geometry for 4/5-mers
    elif n == 6:
        return dict(alpha=0.65, min_k=2, max_diam=2)     # your hexamer settings
    else:
        return dict(alpha=0.65, min_k=2, max_diam=2)

def dense_enough(S: nx.Graph, n: int, overrides: dict | None = None) -> bool:
    p = params_for(n)
    if overrides:
        p.update({k: v for k, v in overrides.items() if v is not None})
    e_needed = math.ceil(p["alpha"] * (n * (n - 1) / 2.0))
    if S.number_of_edges() < e_needed:
        return False
    if min(dict(S.degree()).values()) < p["min_k"]:
        return False
    if p["max_diam"] is not None:
        try:
            if nx.diameter(S) > p["max_diam"]:
                return False
        except nx.NetworkXError:
            return False
    return True

# Optional: relax the "rest are monomers" check.
# Set allow_dimers=True to accept frames with a few stray dimers among the leftovers.
def others_ok(G: nx.Graph, target_nodes: set, allow_dimers: bool = True, max_dimers: int = 1) -> bool:
    others = set(G.nodes) - target_nodes
    if not others:
        return True
    comps = [c for c in nx.connected_components(G.subgraph(others))]
    if not allow_dimers:
        return all(len(c) == 1 for c in comps)
    # allow up to `max_dimers` 2-node components; everything else must be monomers
    dimers = sum(1 for c in comps if len(c) == 2)
    return all(len(c) in (1, 2) for c in comps) and dimers <= max_dimers


def sizes_signature(G: nx.Graph):
    return sorted((len(c) for c in nx.connected_components(G)), reverse=True)

# args
p = argparse.ArgumentParser()
p.add_argument("--topol", required=True)
p.add_argument("--traj",  required=True)
p.add_argument("--target", type=int, required=True)
p.add_argument("--start",  type=int, default=None)   # frame start (inclusive)
p.add_argument("--stop",   type=int, default=None)   # frame stop (exclusive)
p.add_argument("--out",    required=True)            # output XTC
p.add_argument("--alpha", type=float, default=None)  # density; None = sensible default
p.add_argument("--min-res-pairs", type=int, default=None)
p.add_argument("--min-degree", type=int, default=None)
p.add_argument("--max-diam", type=int, default=None)
p.add_argument("--allow-dimers", action="store_true")
p.add_argument("--max-dimers", type=int, default=1)
args = p.parse_args()

# set sensible defaults by n (so 4/5-mers aren’t over-filtered)
# if args.min_res_pairs is None:
#     min_res_pairs = 3 if args.target <= 5 else 4
# else:
#     min_res_pairs = args.min_res_pairs
# if args.alpha is None:
#     alpha = 0.60 if args.target <= 5 else 0.65
# else:
#     alpha = args.alpha
min_res_pairs = 3 if (args.min_res_pairs is None and args.target <= 5) else \
                (4 if args.min_res_pairs is None else args.min_res_pairs)

#topol_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_0-1400ns/start.pdb"
#traj_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_0-1400ns/traj_comp_prot_300K.xtc"

#u = mda.Universe(topol_file, traj_file)
u = mda.Universe(args.topol, args.traj)
peps_ca = [seg.atoms.select_atoms('name CA') for seg in u.segments]
N = len(peps_ca)

over = dict(alpha=args.alpha, min_k=args.min_degree, max_diam=args.max_diam)


#with mda.Writer("/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_0-1400ns/monomers.xtc", u.atoms.n_atoms) as W:
with mda.Writer(args.out, u.atoms.n_atoms) as W:
    #for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
    for ts in tqdm(u.trajectory[args.start:args.stop], desc="Frames"):
        G = nx.Graph(); G.add_nodes_from(range(N))

        # build all edges
        for i, j in combinations(range(N), 2):
            Ai, Aj = peps_ca[i], peps_ca[j]
            ii, jj = capped_distance(
                Ai.positions, Aj.positions, cutoff_ca,
                box=u.dimensions, return_distances=True
            )
            if len(ii) >= min_res_pairs:      # number of CA–CA pairs within cutoff
                G.add_edge(i, j)

        comps = [set(c) for c in nx.connected_components(G)]
        cand = [c for c in comps if len(c) == args.target]
        #cand = [c for c in comps if len(c) == target]

        if len(cand) != 1:
            continue

        H_nodes = cand[0]         # the n-node component
        S = G.subgraph(H_nodes)
        if not dense_enough(S, args.target, overrides=over):
            continue
        if not others_ok(G, H_nodes, allow_dimers=args.allow_dimers, max_dimers=args.max_dimers):
            continue

        # final sanity assertions (optional, helps catch logic slips)
        assert len(H_nodes) == args.target
        assert all(G.degree(n) == 0 for n in set(G.nodes) - H_nodes)

        # DEBUG: log what we’re writing
        print(f"[frame {ts.frame}] cand={cand} hexamer={sorted(H_nodes)} "
              f"edges={S.number_of_edges()} degrees={dict(S.degree())}")

        # H_nodes is the set of segment indices in the accepted n-mer
        hex_atoms = u.atoms[[]]
        for i in sorted(H_nodes):
            hex_atoms = hex_atoms + u.segments[i].atoms

        # center n-mer at box center and wrap molecules whole
        old = u.atoms.positions.copy()
        boxvec = u.dimensions[:3]
        u.atoms.translate(boxvec/2.0 - hex_atoms.center_of_mass())
        #u.atoms.translate(boxvec/2.0 - u.atoms.center_of_mass())
        u.atoms.wrap(compound='segments')   # keeps each peptide whole


        # -------- FINAL sanity check (no PBC): must look like n-mer + monomers --------
        G2 = nx.Graph(); G2.add_nodes_from(range(N))
        for i, j in combinations(range(N), 2):
            Ai, Aj = peps_ca[i], peps_ca[j]
            ii2, jj2 = capped_distance(
                Ai.positions, Aj.positions, cutoff_ca,
                box=None, return_distances=True          # <- Euclidean distances after reimage
            )
            if len(ii2) >= min_res_pairs:
                G2.add_edge(i, j)

        comps = [set(c) for c in nx.connected_components(G2)]
        cand = [c for c in comps if len(c) == args.target]
        if len(cand) != 1:
            u.atoms.positions = old
            continue

        S2 = G2.subgraph(max((set(c) for c in nx.connected_components(G2)), key=len))
        if not dense_enough(S2, args.target, overrides=over):
            u.atoms.positions = old
            continue
        if not others_ok(G2, H_nodes, allow_dimers=args.allow_dimers, max_dimers=args.max_dimers):
            u.atoms.positions = old
            continue

        W.write(u.atoms)
        #W.write(hex_atoms)
        u.atoms.positions = old             # restore for next frame