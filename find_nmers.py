# import MDAnalysis as mda
# from MDAnalysis.lib.distances import capped_distance
# import numpy as np
# import networkx as nx
# from tqdm import tqdm

# target = 6
# cutoff_ca = 8.5
# min_res_pairs = 4

# def dense_enough(S: nx.Graph) -> bool:
#     degs = dict(S.degree())
#     return (min(degs.values()) >= 2) and (S.number_of_edges() >= 12)

# def sizes_signature(G: nx.Graph):
#     return sorted((len(c) for c in nx.connected_components(G)), reverse=True)

# topol_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_1100ns_processed.pdb"
# traj_file  = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_protein_1100-1400ns_300K.xtc"

# u = mda.Universe(topol_file, traj_file)

# # Build an "all CA" AtomGroup with a known peptide ordering
# peps_ca = [seg.atoms.select_atoms('name CA') for seg in u.segments]
# N = len(peps_ca)
# ca_inds_by_pep = [ag.indices for ag in peps_ca]
# all_ca_inds = np.concatenate(ca_inds_by_pep)
# all_ca = u.atoms[all_ca_inds]

# # Map each CA index (in all_ca) to its peptide id
# pep_of_ca = np.empty(all_ca.n_atoms, dtype=np.int16)
# offset = 0
# for p, inds in enumerate(ca_inds_by_pep):
#     n = len(inds)
#     pep_of_ca[offset:offset+n] = p
#     offset += n

# def build_graph_from_counts(counts: np.ndarray, thresh: int) -> nx.Graph:
#     # counts[i,j] = # of CA–CA pairs within cutoff between peptide i and j
#     # (assume upper triangle; symmetrize to be safe)
#     C = counts + counts.T
#     G = nx.Graph()
#     G.add_nodes_from(range(N))
#     ii, jj = np.where(np.triu(C, k=1) >= thresh)
#     G.add_edges_from(zip(ii.tolist(), jj.tolist()))
#     return G

# with mda.Writer("hexamers.xtc", u.atoms.n_atoms) as W:
#     for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
#         # --- One neighbor search across all CA (PBC-aware) ---
#         a, b = capped_distance(
#             all_ca.positions, all_ca.positions, cutoff_ca,
#             box=u.dimensions, return_distances=True
#         )
#         # keep unique pairs and exclude same-atom/self
#         mask = a < b
#         a = a[mask]; b = b[mask]

#         # translate atom pairs -> peptide pairs, discard intra-peptide
#         pi = pep_of_ca[a]; pj = pep_of_ca[b]
#         mask_inter = (pi != pj)
#         pi = pi[mask_inter]; pj = pj[mask_inter]

#         # accumulate contact counts by peptide-pair
#         counts = np.zeros((N, N), dtype=np.int16)
#         np.add.at(counts, (pi, pj), 1)

#         # build graph, enforce 1× hexamer + monomers and density
#         G = build_graph_from_counts(counts, min_res_pairs)
#         sizes = sizes_signature(G)
#         if sizes != [target] + [1]*(N - target):
#             continue

#         H_nodes = max((set(c) for c in nx.connected_components(G)), key=len)
#         S = G.subgraph(H_nodes)
#         if not dense_enough(S):
#             continue

#         # --- Reimage for viewing (same as before, but faster AtomGroup union) ---
#         old = u.atoms.positions.copy()
#         boxvec = u.dimensions[:3]
#         hex_idx = np.concatenate([u.segments[i].atoms.indices for i in sorted(H_nodes)])
#         hex_atoms = u.atoms[hex_idx]
#         u.atoms.translate(boxvec/2.0 - hex_atoms.center_of_mass())
#         u.atoms.wrap(compound='segments')

#         # --- Final no-PBC check using the same all-pairs trick ---
#         a2, b2 = capped_distance(
#             all_ca.positions, all_ca.positions, cutoff_ca,
#             box=None, return_distances=True
#         )
#         m = a2 < b2
#         a2 = a2[m]; b2 = b2[m]
#         pi2 = pep_of_ca[a2]; pj2 = pep_of_ca[b2]
#         inter = (pi2 != pj2); pi2 = pi2[inter]; pj2 = pj2[inter]
#         counts2 = np.zeros((N, N), dtype=np.int16)
#         np.add.at(counts2, (pi2, pj2), 1)
#         G2 = build_graph_from_counts(counts2, min_res_pairs)
#         sizes2 = sizes_signature(G2)
#         if sizes2 != [target] + [1]*(N - target):
#             u.atoms.positions = old; continue
#         S2 = G2.subgraph(max((set(c) for c in nx.connected_components(G2)), key=len))
#         if not dense_enough(S2):
#             u.atoms.positions = old; continue

#         W.write(u.atoms)
#         u.atoms.positions = old


# import MDAnalysis as mda
# from MDAnalysis.lib.distances import capped_distance
# from MDAnalysis.transformations import unwrap, center_in_box, wrap
# from MDAnalysis.topology.guessers import guess_bonds, guess_atom_element
# import networkx as nx
# from itertools import combinations
# from tqdm import tqdm

# topol_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_1100ns.pdb"
# traj_file  = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_protein_1100-1400ns_300K.xtc"
# u = mda.Universe(topol_file, traj_file)
# # Ensure elements exist (improves bond guessing)
# if not hasattr(u.atoms, "elements") or any(e is None for e in getattr(u.atoms, "elements", [])):
#     elems = [guess_atom_element(n) for n in u.atoms.names]
#     u.add_TopologyAttr("elements", elems)
# # Guess covalent bonds (uses vdW radii; respects PBC if box is given)
# pairs = guess_bonds(u.atoms, u.atoms.positions, box=u.dimensions)  # returns list of (i, j) index pairs
# u.add_TopologyAttr("bonds", pairs)
# print(f"Added {len(pairs)} bonds to Universe.")
# peps_ca = [seg.atoms.select_atoms('name CA') for seg in u.segments]
# N = len(peps_ca)

# target = 6
# cutoff_ca = 8.5
# min_res_pairs = 4

# def dense_enough(S: nx.Graph) -> bool:
#     degs = dict(S.degree())
#     return (min(degs.values()) >= 2) and (S.number_of_edges() >= 12)

# def build_graph(peps, cutoff, min_pairs, box=None):
#     """Cα residue-level graph; pass box for PBC-aware, box=None for Euclidean."""
#     G = nx.Graph(); G.add_nodes_from(range(len(peps)))
#     for i, j in combinations(range(len(peps)), 2):
#         Ai, Aj = peps[i], peps[j]
#         ii, jj = capped_distance(Ai.positions, Aj.positions, cutoff,
#                                  box=box, return_distances=True)
#         if len(ii) >= min_pairs:
#             G.add_edge(i, j)
#     return G

# def sizes_signature(G):
#     return sorted((len(c) for c in nx.connected_components(G)), reverse=True)

# def write_centered_hexamer(u, H_nodes, W):
#     # collect atoms in the accepted hexamer
#     hex_atoms = u.atoms[[]]
#     for s in sorted(H_nodes):
#         hex_atoms = hex_atoms + u.segments[s].atoms

#     old = u.atoms.positions.copy()

#     # 1) make coordinates whole (no kwargs; uses bonds/frag info internally)
#     try:
#         _ = unwrap(u.atoms)(u.trajectory.ts)
#     except TypeError:
#         # older/newer variants—if unwrap is unavailable, just skip; next steps still help
#         pass

#     # 2) center box on hexamer COM
#     shift = u.dimensions[:3] / 2.0 - hex_atoms.center_of_mass()
#     u.atoms.translate(shift)

#     # 3) wrap everything back into the unit cell, keeping peptides whole
#     try:
#         u.atoms.wrap(compound='segments')      # keeps each segment together
#     except TypeError:
#         # fallback if 'segments' not supported in your build
#         u.atoms.wrap(compound='fragments')     # needs bonds; otherwise use 'residues'

#     W.write(u.atoms)
#     u.atoms.positions = old

# with mda.Writer("hexamers.xtc", u.atoms.n_atoms) as W:
#     for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
#         # 1) PBC-aware detection on raw coordinates
#         G = build_graph(peps_ca, cutoff_ca, min_res_pairs, box=u.dimensions)
#         sig = sizes_signature(G)
#         print(sig)
#         if sig != [target] + [1]*(N - target):
#             continue

#         H_nodes = max((set(c) for c in nx.connected_components(G)), key=len)
#         S = G.subgraph(H_nodes)
#         if not dense_enough(S):
#             continue

#         # 2) Reimage only the accepted frame for visualization/output
#         write_centered_hexamer(u, H_nodes, W)

#         # 3) Final automated sanity check (no PBC): if it doesn’t look like 6+monomers, drop it
#         #    (Rebuild graph on the reimaged coordinates, Euclidean distances)
#         old = u.atoms.positions.copy()
#         # temporarily put the written positions back to re-check
#         ts = unwrap(u.atoms, compound='segments')(u.trajectory.ts)  # keep whole for check
#         G_check = build_graph(peps_ca, cutoff_ca, min_res_pairs, box=None)
#         ok = (sizes_signature(G_check) == [target] + [1]*(N - target) and
#               dense_enough(G_check.subgraph(H_nodes)))
#         u.atoms.positions = old
#         if not ok:
#             # overwrite last frame by skipping it: simplest is to log and continue;
#             # if you need a strictly clean XTC, write to a list and concat later.
#             print(f"Skipping frame {ts.frame}: failed final no-PBC check.")
#             continue



import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
from itertools import combinations
import networkx as nx
from tqdm import tqdm

target = 6
cutoff_ca = 8.5
min_res_pairs = 4  # your setting

def dense_enough(S: nx.Graph) -> bool:
    degs = dict(S.degree())
    return (min(degs.values()) >= 2) and (S.number_of_edges() >= 12)  # <-- note ()

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
        if not dense_enough(S):
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
        if not dense_enough(S2):
            u.atoms.positions = old
            continue

        W.write(u.atoms)
        u.atoms.positions = old             # restore for next frame

#         # # Optional: reimage so the hexamer is whole and centered (avoids PBC illusions)
#         # old = u.atoms.positions.copy()
#         # box = u.dimensions[:3]
#         # # center on the hexamer COM
#         # hex_atoms = u.atoms[[]]
#         # for i in sorted(H_nodes):
#         #     hex_atoms = hex_atoms + u.segments[i].atoms
#         # shift = box/2.0 - hex_atoms.center_of_mass()
#         # u.atoms.translate(shift)
#         # # wrap molecules back into the box (by segment so each peptide is whole)
#         # u.atoms.wrap(compound='segments')
#         # W.write(u.atoms)
#         # # restore original coordinates for subsequent logic
#         # u.atoms.positions = old


# import MDAnalysis as mda
# from MDAnalysis.lib.distances import capped_distance
# from itertools import combinations
# import networkx as nx
# from tqdm import tqdm

# topol_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_1100ns_processed.pdb"
# traj_file  = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns//traj_comp_prot_pbc_1100-1150ns.xtc"

# u = mda.Universe(topol_file, traj_file)
# peps_ca = [seg.atoms.select_atoms('name CA') for seg in u.segments]
# N = len(peps_ca)

# target = 6
# cutoff_ca = 8.5
# min_res_pairs = 4  # your current setting

# def dense_enough(S: nx.Graph) -> bool:
#     degs = dict(S.degree())
#     return (min(degs.values()) >= 2) and (S.number_of_edges() >= 12)  # your current setting

# with mda.Writer("hexamers.xtc", u.atoms.n_atoms) as W:
#     for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
#         G = nx.Graph(); G.add_nodes_from(range(N))

#         # build all edges first
#         for i, j in combinations(range(N), 2):
#             Ai, Aj = peps_ca[i], peps_ca[j]
#             ii, jj = capped_distance(
#                 Ai.positions, Aj.positions, cutoff_ca,
#                 box=u.dimensions, return_distances=True
#             )
#             if len(ii) >= min_res_pairs:   # number of CA–CA pairs within cutoff
#                 G.add_edge(i, j)

#         # strict target enforcement
#         comps = [set(c) for c in nx.connected_components(G)]
#         sizes = sorted([len(c) for c in comps], reverse=True)

#         # require exactly one cluster of size 'target' and the rest monomers
#         if sizes != [target] + [1]*(N - target):
#             continue

#         H_nodes = max(comps, key=len)        # the unique non-monomer component
#         S = G.subgraph(H_nodes)

#         if not dense_enough(S):
#             continue

#         # final sanity assertions (optional, helps catch logic slips)
#         assert len(H_nodes) == target
#         assert all(G.degree(n) == 0 for n in set(G.nodes) - H_nodes)

#         W.write(u.atoms)


# import MDAnalysis as mda
# from MDAnalysis.lib.distances import capped_distance
# from itertools import combinations
# import networkx as nx
# from tqdm import tqdm

# # Load your system
# topol_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_1100ns_processed.pdb"
# traj_file  = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_pbc_1100-1150ns.xtc"

# u = mda.Universe(topol_file, traj_file)

# # One CA per residue per segment
# peps_ca = [seg.atoms.select_atoms('name CA') for seg in u.segments]
# N = len(peps_ca)

# target = 6
# cutoff_ca = 8.5       # Å
# min_res_pairs = 4     # 3–4 for 11-mers

# def dense_enough(S: nx.Graph) -> bool:
#     degs = dict(S.degree())
#     return (min(degs.values()) >= 2) and (S.number_of_edges() >= 12)

# with mda.Writer("hexamers.xtc", u.atoms.n_atoms) as W:
#     for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
#         G = nx.Graph(); G.add_nodes_from(range(N))

#         # Build edges
#         for i, j in combinations(range(N), 2):
#             Ai, Aj = peps_ca[i], peps_ca[j]
#             ii, jj = capped_distance(
#                 Ai.positions, Aj.positions, cutoff_ca,
#                 box=u.dimensions, return_distances=True
#             )
#             if len(ii) >= min_res_pairs:
#                 G.add_edge(i, j)

#         # Find hexamers and enforce compactness + others are monomers
#         comps = [set(c) for c in nx.connected_components(G)]
#         hexas = [c for c in comps if len(c) == target]

#         if len(hexas) == 1:
#             H_nodes = next(iter(hexas))
#             H = G.subgraph(H_nodes)
#             others = set(G.nodes) - H_nodes
#             if all(G.degree(n) == 0 for n in others) and dense_enough(H):
#                 W.write(u.atoms)


# import MDAnalysis as mda
# from MDAnalysis.analysis.distances import distance_array
# from MDAnalysis.lib.distances import capped_distance
# from itertools import combinations
# import networkx as nx
# import matplotlib.pyplot as plt
# from tqdm import tqdm

# # Load your system
# topol_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_1100ns_processed.pdb"
# traj_file = "/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/traj_comp_prot_pbc_1100-1150ns.xtc"

# u = mda.Universe(topol_file, traj_file)
# peps_ca = [seg.atoms.select_atoms('name CA') for seg in u.segments]
# N = len(peps_ca)
# target = 6

# cutoff_ca = 8.5           # Å
# min_res_pairs = 3         # 3–4 for 11-mers is good

# def dense_enough(S):
#     degs = dict(S.degree())
#     return (min(degs.values()) >= 2) and (S.number_of_edges() >= 10)

# with mda.Writer("hexamers.xtc", u.atoms.n_atoms) as W:
#     for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
#         G = nx.Graph(); G.add_nodes_from(range(N))
#         for i, j in combinations(range(N), 2):
#             Ai, Aj = peps_ca[i], peps_ca[j]
#             ii, jj, _ = capped_distance(Ai.positions, Aj.positions, cutoff_ca,
#                                         box=u.dimensions, return_distances=True)
#             # how many distinct residue pairs are close?
#             if len(set(zip(Ai.resids[ii], Aj.resids[jj]))) >= min_res_pairs:
#                 G.add_edge(i, j)
#             comps = [set(c) for c in nx.connected_components(G)]
#             hexas = [c for c in comps if len(c) == target]

#             if len(hexas) == 1:
#                 H = G.subgraph(next(iter(hexas)))
#                 # all others must be monomers (no edges) and hexamer must be dense
#                 others = set(G.nodes) - set(H.nodes)
#                 if all(G.degree(n) == 0 for n in others) and dense_enough(H):
#                     W.write(u.atoms)




# # pick heavy atoms for each peptide
# peps_heavy = [u.select_atoms(f"segid {seg.segid} and (not name H*)")
#               for seg in u.segments]
# N = len(peps_heavy)
# cutoff = 4.5      # Å for heavy-atom contacts
# min_pairs = 25    # require at least this many atom-atom contacts to form an edge
# target = 6

# with mda.Writer("hexamers.xtc", u.atoms.n_atoms) as W:
#     for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
#         G = nx.Graph(); G.add_nodes_from(range(N))
#         for i, j in combinations(range(N), 2):
#             Ai = peps_heavy[i].positions
#             Aj = peps_heavy[j].positions
#             # returns (idx_i, idx_j, dists); ask only for pairs → faster
#             ii, jj = capped_distance(Ai, Aj, cutoff, box=u.dimensions, return_distances=True)
#             if len(ii) >= min_pairs:
#                 G.add_edge(i, j)

#         comps = [set(c) for c in nx.connected_components(G)]
#         # candidate hexamers with a compactness check (min degree ≥ 2)
#         hexamers = []
#         for c in comps:
#             if len(c) == target:
#                 S = G.subgraph(c)
#                 if min(dict(S.degree()).values(), default=0) >= 2:
#                     hexamers.append(c)

#         # exactly one hexamer; everyone else must be monomers
#         if len(hexamers) == 1:
#             hex_nodes = next(iter(hexamers))
#             others = set(G.nodes) - hex_nodes
#             if all(G.degree(n) == 0 for n in others):
#                 W.write(u.atoms)
#                 # (optional) draw if you want



# peps = [seg.atoms for seg in u.segments]
# print(peps)
# N = len(peps)                 # total # of peptides
# cutoff = 2.5                        # A; adjust to “contact” cutoff
# target = 6                         # e.g. n=4 if you want tetramers
# cnt = 0
# # Prepare an output writer for the selected frames
# with mda.Writer("/home/tank/jsheppard/alpha_syn/NACcore/octamer/Analyses_1100-1400ns/hexamers.xtc", u.atoms.n_atoms) as W:

#     # Loop over frames
#     for ts in tqdm(u.trajectory, desc="Frames", total=u.trajectory.n_frames):
#         # build a graph whose nodes are peptide indices [0…N-1]
#         G = nx.Graph()
#         G.add_nodes_from(range(N))

#         # For each pair of peptides, test if ANY atom–atom contact < cutoff
#         for i, j in combinations(range(N), 2):
#             ci = peps[i].positions     # ∼(natoms_i, 3)
#             cj = peps[j].positions
#             dists = distance_array(ci, cj, box=u.dimensions)    # shape (ni, nj)
#             contacts = (dists < cutoff)                         # Boolean array
#             if contacts.any():
#                 G.add_edge(i, j)

#         # Find connected components → cluster sizes
#         sizes = [len(c) for c in nx.connected_components(G)]

#         # Check “exactly one cluster of size target, rest monomers”
#         if sizes.count(target) == 1 and all(s in (1, target) for s in sizes):
#             cnt += 1
#             if cnt > 50:
#                 pos = nx.spring_layout(G)    # “force‐directed” layout
#                 plt.figure(figsize=(6,6))
#                 nx.draw(
#                     G,
#                     pos,
#                     with_labels=True,
#                     node_size=300,
#                     node_color="skyblue",
#                     edge_color="gray",
#                     font_size=10
#                 )
#                 plt.title(f"Peptide–contact graph, frame {ts.frame}")
#                 plt.axis("off")
#                 plt.show()
#             W.write(u.atoms)   # save this frame
