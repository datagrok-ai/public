#name: PharmacophoreStage3IsolatePocket
#description: Stage 3 - isolate binding-pocket atoms in the consensus frame. Cutoff (default) keeps protein atoms within `pocket_radius` of any ligand atom; DBSCAN (advanced) clusters the union point cloud and keeps the largest non-noise cluster.
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, numpy, pandas, scipy, scikit-learn]
#meta.cache: none
#input: dataframe aligned_structures {caption: aligned_structures from Stage 2a}
#input: string pocket_method = "cutoff" {caption: 'cutoff' (default) or 'dbscan'}
#input: double pocket_radius = 5.0 {caption: Distance cutoff in A (cutoff path only)}
#input: double dbscan_eps = 4.0 {caption: DBSCAN eps in A (DBSCAN path only)}
#input: int dbscan_min_samples = 6 {caption: DBSCAN min_samples (DBSCAN path only)}
#output: dataframe pocket_atoms

# Pure-numpy PDB parsing + scipy.cKDTree cutoff query + sklearn.cluster.DBSCAN.
# Cofactors (HOH, SO4, HEM, NAG, FAD, ...) are filtered before pocket search so they
# can't act as the ligand seed.

import json

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN


# In-script deny-list (kept in sync with files/cofactor-denylist.csv).
# Per blueprint OQ1: the CSV is the documentation source, this set is what runs.
COFACTOR_DENYLIST = {
    'HOH', 'DOD',
    'SO4', 'PO4',
    'CL', 'NA', 'K', 'CA', 'MG', 'ZN', 'FE', 'MN', 'NI', 'CU',
    'EDO', 'GOL', 'PEG', 'PG4', 'PGE',
    'DMS', 'ACT', 'TRS', 'MES', 'IMD',
    'HEM', 'HEC',
    'NAG', 'BMA', 'MAN', 'FUC', 'GAL',
    'FAD', 'NAD', 'NAP', 'ADP', 'ATP', 'GDP', 'GTP', 'COA', 'FMN', 'SAM', 'SAH',
}


def parse_pdb_atoms(pdb_text):
    """
    Parse ATOM/HETATM lines from a PDB block. Returns a list of dicts:
    {record, chain, res_name, res_seq, atom_name, element, x, y, z}.
    """
    out = []
    for line in pdb_text.split('\n'):
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        if len(line) < 54:
            continue
        try:
            chain = line[21]
            res_name = line[17:20].strip()
            res_seq = int(line[22:26].strip())
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) >= 78 else atom_name[:1]
        except ValueError:
            continue
        out.append({
            'record': 'ATOM' if line.startswith('ATOM') else 'HETATM',
            'chain': chain, 'res_name': res_name, 'res_seq': res_seq,
            'atom_name': atom_name, 'element': element or atom_name[:1],
            'x': x, 'y': y, 'z': z,
        })
    return out


def apply_transform(coords, transform_flat):
    """Apply a 16-element row-major 4x4 transform to an Nx3 coord array."""
    T = np.array(transform_flat, dtype=float).reshape(4, 4)
    R = T[:3, :3]
    t = T[:3, 3]
    return (coords @ R.T) + t


def atoms_within_radius_of_any_ligand(protein_coords, ligand_coords, radius):
    """
    Return a boolean array (len = len(protein_coords)) marking atoms within `radius`
    of any ligand atom. Built on cKDTree.query for O(N log M) instead of N*M.
    """
    if ligand_coords.shape[0] == 0 or protein_coords.shape[0] == 0:
        return np.zeros(protein_coords.shape[0], dtype=bool)
    tree = cKDTree(ligand_coords)
    distances, _ = tree.query(protein_coords, k=1, distance_upper_bound=radius)
    return np.isfinite(distances)


# --- Main body -----------------------------------------------------------

if not isinstance(aligned_structures, pd.DataFrame) or len(aligned_structures) == 0:
    raise ValueError('Stage 3: aligned_structures is empty.')

method = (pocket_method or 'cutoff').strip().lower()
if method not in ('cutoff', 'dbscan'):
    print(f'Stage 3: unknown pocket_method "{pocket_method}"; defaulting to cutoff.')
    method = 'cutoff'

rows = []
for _, row in aligned_structures.iterrows():
    pdb_id = str(row['pdb_id'])
    pdb_text = row.get('original_pdb') or ''
    transform_json = row.get('transform_4x4_json') or '[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]'
    transform = json.loads(transform_json)

    atoms = parse_pdb_atoms(pdb_text)
    if not atoms:
        print(f'Stage 3: {pdb_id} has no parseable atoms')
        continue

    raw_coords = np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=float)
    consensus_coords = apply_transform(raw_coords, transform)

    ligand_idx, protein_idx = [], []
    for i, a in enumerate(atoms):
        if a['record'] == 'HETATM':
            if a['res_name'].upper() in COFACTOR_DENYLIST:
                continue
            ligand_idx.append(i)
        else:
            protein_idx.append(i)

    if not ligand_idx:
        print(f'Stage 3: {pdb_id} has no non-cofactor ligand HETATMs — skipping')
        continue
    if not protein_idx:
        print(f'Stage 3: {pdb_id} has no protein ATOM records — skipping')
        continue

    ligand_xyz = consensus_coords[ligand_idx]
    protein_xyz = consensus_coords[protein_idx]
    within = atoms_within_radius_of_any_ligand(protein_xyz, ligand_xyz, pocket_radius)

    n_pocket = int(within.sum())
    print(f'Stage 3: {pdb_id} pocket={n_pocket} (of {len(protein_idx)} protein atoms); '
          f'{len(ligand_idx)} ligand atoms')

    # A residue is "in the pocket" if any of its atoms is within `pocket_radius`
    # of any ligand atom. We emit only the Calpha of those residues — that's all
    # Stage 2b (pocket-only Kabsch) and the Mol* overlay need; side chains were
    # dead weight in the original `all atoms` output.
    pocket_residues = set()
    for j, mask_hit in enumerate(within):
        if not mask_hit:
            continue
        a = atoms[protein_idx[j]]
        pocket_residues.add((a['chain'], a['res_seq']))

    n_pocket_residues = len(pocket_residues)
    print(f'Stage 3: {pdb_id} pocket spans {n_pocket_residues} residue(s) '
          f'(emitting only Calpha + ligand atoms)')

    for j in range(len(protein_idx)):
        i = protein_idx[j]
        a = atoms[i]
        if (a['chain'], a['res_seq']) not in pocket_residues:
            continue
        if a['atom_name'] != 'CA':
            continue
        rows.append({
            'pdb_id': pdb_id, 'chain': a['chain'], 'res_name': a['res_name'],
            'res_seq': a['res_seq'], 'atom_name': a['atom_name'],
            'element': a['element'],
            'x_consensus': float(consensus_coords[i, 0]),
            'y_consensus': float(consensus_coords[i, 1]),
            'z_consensus': float(consensus_coords[i, 2]),
            'ligand_seed': False,
            'cluster_label': 0,
        })

    # Always include the ligand atoms themselves (so downstream stages have them)
    for i in ligand_idx:
        a = atoms[i]
        rows.append({
            'pdb_id': pdb_id, 'chain': a['chain'], 'res_name': a['res_name'],
            'res_seq': a['res_seq'], 'atom_name': a['atom_name'],
            'element': a['element'],
            'x_consensus': float(consensus_coords[i, 0]),
            'y_consensus': float(consensus_coords[i, 1]),
            'z_consensus': float(consensus_coords[i, 2]),
            'ligand_seed': True,
            'cluster_label': 0,
        })

# Optional DBSCAN re-cluster over the full pocket cloud (consensus frame).
# We cluster only the pocket protein atoms (ligand seeds stay in regardless).
if method == 'dbscan' and rows:
    seed_rows = [r for r in rows if r['ligand_seed']]
    cand_rows = [r for r in rows if not r['ligand_seed']]
    if cand_rows:
        cloud = np.array([[r['x_consensus'], r['y_consensus'], r['z_consensus']]
                          for r in cand_rows], dtype=float)
        db = DBSCAN(eps=float(dbscan_eps), min_samples=int(dbscan_min_samples))
        labels = db.fit_predict(cloud)
        non_noise = labels[labels >= 0]
        if non_noise.size > 0:
            unique, counts = np.unique(non_noise, return_counts=True)
            winner = int(unique[int(np.argmax(counts))])
            kept = []
            for r, lbl in zip(cand_rows, labels):
                if int(lbl) == winner:
                    r['cluster_label'] = winner
                    kept.append(r)
            print(f'Stage 3: DBSCAN(eps={dbscan_eps}, min_samples={dbscan_min_samples}) '
                  f'kept cluster {winner} with {len(kept)} of {len(cand_rows)} atoms; '
                  f'noise dropped.')
            rows = kept + seed_rows
        else:
            print('Stage 3: DBSCAN labeled every point as noise — falling back to cutoff result.')

pocket_atoms = pd.DataFrame(rows)
print(f'Stage 3: emitting {len(pocket_atoms)} pocket atoms ({method} path)')
