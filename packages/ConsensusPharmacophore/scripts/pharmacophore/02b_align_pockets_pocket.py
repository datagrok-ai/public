#name: PharmacophoreStage2bAlignPocketsPocket
#description: Stage 2b - pocket-only Calpha Kabsch refinement. Replaces pass-1 transforms with pass-2 over the consensus-pocket Calpha subset. Filters by selected cluster (selected_pdb_ids).
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, numpy, pandas]
#meta.cache: none
#input: dataframe aligned_structures {caption: Stage 2a output (with pass-1 transforms)}
#input: dataframe pocket_atoms {caption: Stage 3 output (pocket Calpha atoms + ligand seeds)}
#input: string selected_pdb_ids = "[]" {caption: JSON array of PDB IDs to keep}
#output: dataframe aligned_structures_v2 {caption: Pass-2 transforms over pocket Calpha}

import json
import math
import numpy as np
import pandas as pd

IDENTITY_4X4 = np.eye(4).flatten().tolist()
_MIN_COMMON_CA = 5


def parse_ca_by_resid(pdb_text):
    """{ (chain, resid) -> (x, y, z) } for ATOM Calpha records only."""
    out = {}
    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM') or len(line) < 54:
            continue
        if line[12:16].strip() != 'CA':
            continue
        try:
            chain = line[21]
            resid = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        out.setdefault((chain, resid), (x, y, z))
    return out


def kabsch(P, Q):
    """Return R such that R @ P_col ~ Q_col (P, Q centered Nx3 arrays)."""
    H = P.T @ Q
    U, _S, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def apply_transform_to_points(points, T):
    """Apply 4x4 row-major transform to Nx3 array. Returns Nx3."""
    R = T[:3, :3]
    t = T[:3, 3]
    return (points @ R.T) + t


# --- Main body -----------------------------------------------------------

if not isinstance(aligned_structures, pd.DataFrame) or len(aligned_structures) == 0:
    raise ValueError('Stage 2b: aligned_structures is empty.')
if not isinstance(pocket_atoms, pd.DataFrame) or len(pocket_atoms) == 0:
    raise ValueError('Stage 2b: pocket_atoms is empty.')

# Parse the cluster picker's selected PDB list.
try:
    selected = set(json.loads(selected_pdb_ids or '[]'))
except (json.JSONDecodeError, TypeError):
    selected = set()
use_all = (not selected) or ('all' in selected)

# Pocket residue set per (pdb_id, chain, resid). Exclude ligand-seed atoms — only
# protein Calpha residues define the pocket alignment subset.
pocket_keys = set()
seed_avail = 'ligand_seed' in pocket_atoms.columns
for i in range(len(pocket_atoms)):
    row = pocket_atoms.iloc[i]
    if seed_avail and bool(row['ligand_seed']):
        continue
    if str(row.get('atom_name', '')).strip() != 'CA':
        continue
    pocket_keys.add((str(row['pdb_id']), str(row['chain']), int(row['res_seq'])))

print(f'Stage 2b: {len(pocket_keys)} unique pocket-Calpha keys across all PDBs')

# Template = PDB with smallest pass-1 rmsd_to_ref (usually 0, the template).
template_idx = aligned_structures['rmsd_to_ref'].astype(float).idxmin()
template_id = str(aligned_structures.at[template_idx, 'pdb_id'])
template_text = str(aligned_structures.at[template_idx, 'original_pdb'])
template_ca_all = parse_ca_by_resid(template_text)

# Template's pocket Calpha set (consensus-pocket residues defined by the template).
template_pocket_keys = {(c, r) for (pid, c, r) in pocket_keys if pid == template_id}
template_pocket_ca = {k: template_ca_all[k] for k in template_pocket_keys if k in template_ca_all}
print(f'Stage 2b: template {template_id} contributes {len(template_pocket_ca)} pocket Calpha atoms')

out_rows = []
bootstrap_failed_count = 0

for i in range(len(aligned_structures)):
    src = aligned_structures.iloc[i].to_dict()
    pid = str(src['pdb_id'])

    if not use_all and pid not in selected:
        continue

    if pid == template_id:
        src['transform_4x4_json'] = json.dumps(IDENTITY_4X4)
        src['rmsd_to_ref'] = 0.0
        src['bootstrap_failed'] = False
        out_rows.append(src)
        continue

    pdb_text = str(src.get('original_pdb', '') or '')
    if not pdb_text:
        src['bootstrap_failed'] = True
        out_rows.append(src)
        bootstrap_failed_count += 1
        continue

    # Mobile pocket Calpha set: residues that ARE in this PDB's pocket atoms AND have a Calpha record.
    mobile_pocket_keys = {(c, r) for (mpid, c, r) in pocket_keys if mpid == pid}
    mobile_ca_all = parse_ca_by_resid(pdb_text)
    mobile_pocket_ca = {k: mobile_ca_all[k] for k in mobile_pocket_keys if k in mobile_ca_all}

    # Common keys = residues in BOTH mobile's and template's pocket. Sort for stable order.
    common = sorted(set(mobile_pocket_ca.keys()) & set(template_pocket_ca.keys()))
    if len(common) < _MIN_COMMON_CA:
        print(f'Stage 2b: {pid} has only {len(common)} common pocket Calpha — bootstrap_failed')
        src['bootstrap_failed'] = True
        out_rows.append(src)
        bootstrap_failed_count += 1
        continue

    P = np.array([mobile_pocket_ca[k] for k in common], dtype=float)
    Q = np.array([template_pocket_ca[k] for k in common], dtype=float)

    # Pass-1 RMSD over the SAME pocket subset (apply pass-1 transform to mobile, measure vs template).
    try:
        T1 = np.array(json.loads(src['transform_4x4_json']), dtype=float).reshape(4, 4)
        P_pass1 = apply_transform_to_points(P, T1)
        pass1_pocket_rmsd = float(np.sqrt(np.mean(np.sum((P_pass1 - Q) ** 2, axis=1))))
    except (ValueError, KeyError, json.JSONDecodeError):
        pass1_pocket_rmsd = float('nan')

    # Pass-2 Kabsch from raw mobile -> template (template is in consensus = original frame).
    pc = P.mean(axis=0)
    qc = Q.mean(axis=0)
    R = kabsch(P - pc, Q - qc)
    rotated = (R @ (P - pc).T).T
    pass2_pocket_rmsd = float(np.sqrt(np.mean(np.sum((rotated - (Q - qc)) ** 2, axis=1))))
    t = qc - R @ pc

    T2 = np.eye(4)
    T2[:3, :3] = R
    T2[:3, 3] = t
    if not np.isfinite(T2).all():
        print(f'Stage 2b: {pid} produced non-finite transform — bootstrap_failed')
        src['bootstrap_failed'] = True
        out_rows.append(src)
        bootstrap_failed_count += 1
        continue

    # By construction pass-2 IS the Kabsch optimum on this subset, so pass-2 RMSD <= pass-1 RMSD.
    # If somehow it's worse (numerical edge case), fall back to pass-1.
    if (math.isfinite(pass1_pocket_rmsd) and pass2_pocket_rmsd > pass1_pocket_rmsd + 1e-6):
        print(f'Stage 2b: {pid} pass-2 RMSD {pass2_pocket_rmsd:.3f} > pass-1 {pass1_pocket_rmsd:.3f} — '
              'keeping pass-1 transform; bootstrap_failed')
        src['bootstrap_failed'] = True
        out_rows.append(src)
        bootstrap_failed_count += 1
        continue

    print(f'Stage 2b: {pid} pass-1-pocket={pass1_pocket_rmsd:.3f}A -> pass-2-pocket={pass2_pocket_rmsd:.3f}A '
          f'(common Calpha: {len(common)})')
    src['transform_4x4_json'] = json.dumps(T2.flatten().tolist())
    src['rmsd_to_ref'] = pass2_pocket_rmsd
    src['bootstrap_failed'] = False
    out_rows.append(src)

aligned_structures_v2 = pd.DataFrame(out_rows)
print(f'Stage 2b: emitted {len(aligned_structures_v2)} rows; '
      f'{bootstrap_failed_count} bootstrap failures')
