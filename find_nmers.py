import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import math
import numpy as np
from itertools import combinations
import networkx as nx
from tqdm import tqdm

target = 6
cutoff_ca = 8.5
min_res_pairs = 4  # your setting

def dense_enough(S: nx.Graph, n: int, alpha: float = 0.65, min_k: int = 2, use_diameter: bool = True) -> bool:
    """
    Accept the component if:
      - it has at least alpha * (n choose 2) edges (density threshold),
      - every node has degree >= min_k (no chains),
      - optional: graph diameter <= 2 (prevents stringy shapes, but allows rings).
    """
    if len(S) != n:
        return False

    # edge density threshold
    e_required = math.ceil(alpha * (n * (n - 1) / 2.0))
    if S.number_of_edges() < e_required:
        return False

    # minimum degree
    if min(dict(S.degree()).values()) < min_k:
        return False

    if use_diameter:
        try:
            if nx.diameter(S) > 2:
                return False
        except nx.NetworkXError:
            # should not happen (component is connected), but be safe
            return False

    return True

def sizes_signature(G: nx.Graph):
    return sorted((len(c) for c in nx.connected_components(G)), reverse=True)

topol_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_1100ns_processed.pdb"
traj_file  = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_protein_1100-1400ns_300K.xtc"

u = mda.Universe(topol_file, traj_file)
peps_ca = [seg.atoms.select_atoms('name CA') for seg in u.segments]
N = len(peps_ca)

with mda.Writer("hexamers.xtc", u.atoms.n_atoms) as W:
    for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
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
        sizes = sorted([len(c) for c in comps], reverse=True)

        # strict requirement: exactly one n-mer + (N-target) monomers
        if sizes != [target] + [1]*(N - target):
            continue

        H_nodes = max(comps, key=len)         # the n-node component
        S = G.subgraph(H_nodes)
        if not dense_enough(S, target):
            continue

        # final sanity assertions (optional, helps catch logic slips)
        assert len(H_nodes) == target
        assert all(G.degree(n) == 0 for n in set(G.nodes) - H_nodes)

        # DEBUG: log what we’re writing
        print(f"[frame {ts.frame}] sizes={sizes} hexamer={sorted(H_nodes)} "
              f"edges={S.number_of_edges()} degrees={dict(S.degree())}")

        # H_nodes is the set of segment indices in the accepted n-mer
        hex_atoms = u.atoms[[]]
        for i in sorted(H_nodes):
            hex_atoms = hex_atoms + u.segments[i].atoms

        # center n-mer at box center and wrap molecules whole
        old = u.atoms.positions.copy()
        boxvec = u.dimensions[:3]
        u.atoms.translate(boxvec/2.0 - hex_atoms.center_of_mass())
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

        sizes2 = sizes_signature(G2)
        if sizes2 != [target] + [1]*(N - target):
            u.atoms.positions = old
            continue

        S2 = G2.subgraph(max((set(c) for c in nx.connected_components(G2)), key=len))
        if not dense_enough(S2, target):
            u.atoms.positions = old
            continue

        W.write(u.atoms)
        u.atoms.positions = old             # restore for next frame