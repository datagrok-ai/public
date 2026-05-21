#name: PharmacophoreStage2AlignPockets
#description: Stage 2a - global Calpha Kabsch alignment of input PDBs. Template is the highest-resolution entry; transform is row-major 4x4 in JSON.
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, numpy, pandas, requests]
#meta.cache: none
#input: dataframe pdb_qc {caption: PDB QC table from Stage 1}
#input: string chain_selection = "auto" {caption: Chain selection ('auto' picks the chain with the most Calpha atoms)}
#input: bool alt_loc_filter = true {caption: Drop altLoc != ' '/'A'}
#input: double min_occupancy = 0.5 {caption: Minimum occupancy}
#output: dataframe aligned_structures

# Pure-numpy implementation — no MDAnalysis dependency.
# PDB v3.3 column layout (1-based, inclusive):
#   7-11   atom serial
#   13-16  atom name      (e.g. ' CA ')
#   17     altLoc
#   18-20  resName
#   22     chainID
#   23-26  resSeq
#   31-38  x  (%8.3f)
#   39-46  y  (%8.3f)
#   47-54  z  (%8.3f)
#   55-60  occupancy

import json
import math

import numpy as np
import pandas as pd
import requests


IDENTITY_4X4 = np.eye(4).flatten().tolist()
_MIN_COMMON_CA = 5


def fetch_pdb(pdb_id):
    """Download a PDB text from RCSB. Raises on HTTP error."""
    resp = requests.get(f'https://files.rcsb.org/download/{pdb_id}.pdb', timeout=30)
    resp.raise_for_status()
    return resp.text


def extract_ca_positions(pdb_text, alt_loc_filter_v=True, min_occupancy_v=0.5):
    """
    Parse a PDB text and return {chainID: {resid: (x, y, z)}} for Calpha atoms,
    after altLoc / occupancy QC. Multi-character residue numbers and standard
    chains are handled per PDB v3.3 fixed-width spec.
    """
    by_chain = {}
    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM'):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name != 'CA':
            continue
        altloc = line[16]
        if alt_loc_filter_v and altloc not in (' ', 'A'):
            continue
        try:
            occ = float(line[54:60].strip()) if len(line) >= 60 else 1.0
        except ValueError:
            occ = 1.0
        if occ < min_occupancy_v:
            continue
        chain_id = line[21]
        try:
            resid = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        ch = by_chain.setdefault(chain_id, {})
        # First-occurrence wins (handles duplicate residues from altLocs that slipped through)
        ch.setdefault(resid, (x, y, z))
    return by_chain


def best_chain(ca_by_chain):
    """Chain with the most Calpha atoms; None if empty."""
    if not ca_by_chain:
        return None
    return max(ca_by_chain.keys(), key=lambda c: len(ca_by_chain[c]))


def kabsch(P, Q):
    """
    Return the 3x3 rotation matrix R such that `R @ P_col ~= Q_col`.
    P and Q are N x 3 arrays of CENTERED points (zero centroid).
    Uses the classic Kabsch SVD with sign correction.
    """
    H = P.T @ Q  # 3 x 3 covariance
    U, _S, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def align(mobile_ca, ref_ca, chain_sel='auto'):
    """
    Compute the 4x4 transform that maps mobile -> ref via shared resids.
    Returns (T_row_major_list, rmsd) or None on failure.
    """
    if chain_sel == 'auto':
        mc = best_chain(mobile_ca)
        rc = best_chain(ref_ca)
    else:
        mc = rc = chain_sel
    if not mc or not rc or mc not in mobile_ca or rc not in ref_ca:
        return None

    mobile_dict = mobile_ca[mc]
    ref_dict = ref_ca[rc]
    common = sorted(set(mobile_dict.keys()) & set(ref_dict.keys()))
    if len(common) < _MIN_COMMON_CA:
        return None

    mobile_pos = np.array([mobile_dict[r] for r in common], dtype=float)
    ref_pos = np.array([ref_dict[r] for r in common], dtype=float)
    mc_centroid = mobile_pos.mean(axis=0)
    rc_centroid = ref_pos.mean(axis=0)
    P = mobile_pos - mc_centroid
    Q = ref_pos - rc_centroid

    R = kabsch(P, Q)
    # rmsd over the matched set after applying R
    rotated = (R @ P.T).T
    rmsd = float(np.sqrt(np.mean(np.sum((rotated - Q) ** 2, axis=1))))

    t = rc_centroid - R @ mc_centroid

    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t
    if not np.isfinite(T).all():
        return None
    return T.flatten().tolist(), rmsd


def _res_key(v):
    try:
        return float(v) if v is not None and not (isinstance(v, float) and math.isnan(v)) else math.inf
    except (TypeError, ValueError):
        return math.inf


# --- Main script body ----------------------------------------------------

if not isinstance(pdb_qc, pd.DataFrame) or len(pdb_qc) == 0:
    raise ValueError('Stage 2a: pdb_qc is empty.')

unique_df = pdb_qc.drop_duplicates(subset='pdb_id', keep='first').reset_index(drop=True).copy()
unique_df['_res_key'] = unique_df['resolution'].apply(_res_key)
template_idx = unique_df['_res_key'].idxmin()
template_id = str(unique_df.at[template_idx, 'pdb_id'])
print(f'Stage 2a: {len(unique_df)} unique PDB(s); template = {template_id}')

# Fetch each PDB once.
pdb_texts = {}
for pid in unique_df['pdb_id'].astype(str):
    try:
        pdb_texts[pid] = fetch_pdb(pid)
        print(f'Stage 2a: fetched {pid} ({len(pdb_texts[pid])} bytes)')
    except Exception as e:
        print(f'Stage 2a: failed to fetch {pid}: {e}')
        pdb_texts[pid] = ''

# Extract Calpha positions for each.
ca_by_pdb = {}
for pid, text in pdb_texts.items():
    if not text:
        continue
    ca_by_pdb[pid] = extract_ca_positions(text, alt_loc_filter, min_occupancy)
    n_ca = sum(len(c) for c in ca_by_pdb[pid].values())
    print(f'Stage 2a: {pid} has {n_ca} Calpha atom(s) across {len(ca_by_pdb[pid])} chain(s)')

ref_ca = ca_by_pdb.get(template_id)
if ref_ca is None or not ref_ca:
    raise ValueError(f'Stage 2a: template {template_id} has no usable Calpha atoms.')

# Align each non-template against the template.
rows = []
for pid in unique_df['pdb_id'].astype(str):
    pdb_text = pdb_texts.get(pid, '')
    if pid == template_id:
        rows.append({
            'pdb_id': pid,
            'original_pdb': pdb_text,
            'transform_4x4_json': json.dumps(IDENTITY_4X4),
            'rmsd_to_ref': 0.0,
            'ref_pdb_id': template_id,
            'bootstrap_failed': False,
            'cluster_id': 'C1',
        })
        continue
    mobile_ca = ca_by_pdb.get(pid)
    res = align(mobile_ca, ref_ca, chain_selection) if mobile_ca else None
    if res is None:
        print(f'Stage 2a: alignment failed for {pid}')
        rows.append({
            'pdb_id': pid, 'original_pdb': pdb_text,
            'transform_4x4_json': json.dumps(IDENTITY_4X4),
            'rmsd_to_ref': float('nan'),
            'ref_pdb_id': template_id, 'bootstrap_failed': True, 'cluster_id': 'C1',
        })
        continue
    T_flat, rmsd = res
    print(f'Stage 2a: {pid} aligned, rmsd={rmsd:.3f} A')
    rows.append({
        'pdb_id': pid, 'original_pdb': pdb_text,
        'transform_4x4_json': json.dumps(T_flat),
        'rmsd_to_ref': rmsd,
        'ref_pdb_id': template_id, 'bootstrap_failed': False, 'cluster_id': 'C1',
    })

aligned_structures = pd.DataFrame(rows)
print(f'Stage 2a: returning {len(aligned_structures)} rows')
