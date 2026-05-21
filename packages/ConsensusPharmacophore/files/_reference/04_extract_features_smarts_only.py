# ARCHIVED — original Stage 4 (SMARTS-only, ligand atoms only).
#
# Production Stage 4 has been replaced by the ProLIF-based interaction
# extractor in scripts/pharmacophore/04_extract_features.py, which emits
# only features that actually contact the protein. Keep this file for:
#   - regression comparison (ligand-only vs interaction-only feature counts)
#   - falling back if ProLIF is unavailable in the target environment
#   - reading the SMARTS-family algorithm in its simplest form
#
# To resurrect: rename back to scripts/pharmacophore/04_extract_features.py
# and add the smarts_csv argument back to the orchestrator's Stage 4 call.
#
# Function name and inputs preserved below for reference; the leading "#name:"
# has been demoted to "# name:" so Datagrok won't register this file.
#
# name: PharmacophoreStage4ExtractFeatures
# description: Stage 4 - per-ligand pharmacophore feature extraction. AssignBondOrdersFromTemplate (with rdFMCS fallback) + SMARTS-family dictionary (D/A/a/H/P/N/X). Centroid each match in 3D and lift to consensus frame.
# language: python
# environment: channels: [conda-forge, defaults], dependencies: [python=3.11, numpy, pandas, requests, rdkit]
# meta.cache: none
# input: dataframe aligned_structures {caption: Stage 2b pass-2 aligned structures}
# input: dataframe pocket_atoms {caption: Stage 3 pocket atoms (used to discover ligand comp_ids per PDB)}
# input: string selected_pdb_ids = "[]" {caption: JSON array of PDB IDs to keep}
# input: string smarts_csv {caption: pharmacophore-features.csv content (4 cols feature_id family feature_name smarts)}
# output: dataframe ligand_features {caption: One row per SMARTS match centroid}

import csv
import io
import json
import math

import numpy as np
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS


# Single-letter family code (Datagrok pharmacophore-features.csv) -> canonical name
# used downstream by Stage 5a / Stage 5b / renderer / family-map.ts. Kept in sync
# with src/family-map.ts.
FAMILY_NAME = {
    'D': 'Donor',
    'A': 'Acceptor',
    'a': 'Aromatic',
    'H': 'Hydrophobic',
    'P': 'Positive',
    'N': 'Negative',
    'X': 'Halogen',
}

_MIN_MCS_COVERAGE = 0.7   # rdFMCS fallback rejects matches covering <70% of template atoms


def parse_smarts_csv(csv_text):
    """Returns list of (family_code, smarts_str). Skips header + invalid rows."""
    out = []
    reader = csv.DictReader(io.StringIO(csv_text))
    for row in reader:
        fam = (row.get('family') or '').strip()
        smarts = (row.get('smarts') or '').strip()
        if fam and smarts and fam in FAMILY_NAME:
            out.append((fam, smarts))
    return out


def extract_ligand_block(pdb_text, comp_id):
    """Return PDB block containing only HETATM records for the given comp_id."""
    cid = comp_id.upper()
    lines = []
    for line in pdb_text.split('\n'):
        if line.startswith('HETATM') and line[17:20].strip().upper() == cid:
            lines.append(line)
    lines.append('END')
    return '\n'.join(lines)


def fetch_ccd_smiles(comp_id):
    """Fetch canonical SMILES from RCSB CCD endpoint; prefer SMILES_stereo over SMILES."""
    url = f'https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}'
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    desc = data.get('rcsb_chem_comp_descriptor') or {}
    smiles = desc.get('SMILES_stereo') or desc.get('SMILES')
    if smiles:
        return smiles
    # Fallback: scan pdbx_chem_comp_descriptor array
    for d in data.get('pdbx_chem_comp_descriptor', []):
        if (d.get('type') == 'SMILES_CANONICAL') and d.get('descriptor'):
            return d['descriptor']
    return None


def correct_bond_orders(template_smiles, mobile_mol):
    """
    Apply RDKit's AssignBondOrdersFromTemplate; on AtomValenceException
    (or similar) fall back to rdFMCS substructure projection. Returns
    (corrected_mol, skip_reason_or_None).
    """
    template = Chem.MolFromSmiles(template_smiles)
    if template is None:
        return None, 'template_smiles_invalid'
    template_no_h = Chem.RemoveHs(template)
    try:
        mobile_no_h = Chem.RemoveHs(mobile_mol, sanitize=False)
    except Exception:
        mobile_no_h = mobile_mol

    try:
        corrected = AllChem.AssignBondOrdersFromTemplate(template_no_h, mobile_no_h)
        Chem.SanitizeMol(corrected)
        return corrected, None
    except (ValueError, RuntimeError, Chem.AtomValenceException) as e:
        # rdFMCS fallback. Match the maximum common substructure between
        # template and mobile, project the template's bond orders onto that
        # subset, and continue. Coverage below _MIN_MCS_COVERAGE is rejected.
        try:
            mcs = rdFMCS.FindMCS([template_no_h, mobile_no_h], timeout=10)
            if mcs.numAtoms <= 0 or mcs.numAtoms < _MIN_MCS_COVERAGE * template_no_h.GetNumAtoms():
                return None, f'mcs_coverage_below_threshold (atoms={mcs.numAtoms})'
            patt = Chem.MolFromSmarts(mcs.smartsString)
            mob_match = mobile_no_h.GetSubstructMatch(patt)
            if not mob_match:
                return None, 'mcs_no_match_in_mobile'
            # Build sub-mobile from the matched atoms only and try again
            sub_mobile = Chem.PathToSubmol(mobile_no_h, list(mob_match))
            try:
                corrected = AllChem.AssignBondOrdersFromTemplate(patt, sub_mobile)
                Chem.SanitizeMol(corrected)
                return corrected, None
            except Exception as e2:
                return None, f'fmcs_fallback_failed: {e2}'
        except Exception as e2:
            return None, f'fmcs_engine_failed: {e2}'


def apply_4x4_to_point(T, x, y, z):
    """T is a 4x4 numpy array; returns (nx, ny, nz)."""
    R = T[:3, :3]
    t = T[:3, 3]
    p = np.array([x, y, z], dtype=float)
    out = R @ p + t
    return float(out[0]), float(out[1]), float(out[2])


# --- Main body -----------------------------------------------------------

if not isinstance(aligned_structures, pd.DataFrame) or len(aligned_structures) == 0:
    raise ValueError('Stage 4: aligned_structures is empty.')
if not isinstance(pocket_atoms, pd.DataFrame) or len(pocket_atoms) == 0:
    raise ValueError('Stage 4: pocket_atoms is empty.')

try:
    selected = set(json.loads(selected_pdb_ids or '[]'))
except (json.JSONDecodeError, TypeError):
    selected = set()
use_all = (not selected) or ('all' in selected)

smarts_rows = parse_smarts_csv(smarts_csv or '')
if not smarts_rows:
    raise ValueError('Stage 4: smarts_csv is empty or unparseable.')
print(f'Stage 4: {len(smarts_rows)} SMARTS patterns loaded')

# Discover each PDB's ligand comp_ids from pocket_atoms (where ligand_seed=True).
ligand_by_pdb = {}
for i in range(len(pocket_atoms)):
    row = pocket_atoms.iloc[i]
    if not bool(row.get('ligand_seed')):
        continue
    pid = str(row['pdb_id'])
    comp = str(row.get('res_name', '')).strip().upper()
    if comp:
        ligand_by_pdb.setdefault(pid, set()).add(comp)

# Cache CCD SMILES across PDBs (the same comp_id often appears multiple times).
smiles_cache = {}
def get_smiles(comp_id):
    if comp_id in smiles_cache:
        return smiles_cache[comp_id]
    try:
        s = fetch_ccd_smiles(comp_id)
    except Exception as e:
        print(f'Stage 4: CCD fetch failed for {comp_id}: {e}')
        s = None
    smiles_cache[comp_id] = s
    return s


out_rows = []
for i in range(len(aligned_structures)):
    src = aligned_structures.iloc[i].to_dict()
    pid = str(src['pdb_id'])
    if not use_all and pid not in selected:
        continue
    if pid not in ligand_by_pdb:
        print(f'Stage 4: {pid} has no ligand in pocket_atoms — skip')
        continue

    pdb_text = str(src.get('original_pdb', '') or '')
    if not pdb_text:
        continue
    try:
        T = np.array(json.loads(src['transform_4x4_json']), dtype=float).reshape(4, 4)
    except (ValueError, KeyError, json.JSONDecodeError):
        T = np.eye(4)

    for comp_id in sorted(ligand_by_pdb[pid]):
        ligand_block = extract_ligand_block(pdb_text, comp_id)
        mobile = Chem.MolFromPDBBlock(ligand_block, sanitize=False, removeHs=False)
        if mobile is None or mobile.GetNumAtoms() == 0:
            print(f'Stage 4: {pid}:{comp_id} HETATM block did not parse')
            out_rows.append({
                'pdb_id': pid, 'ligand_comp_id': comp_id,
                'family': '', 'feature_idx': -1, 'atom_indices_json': '[]',
                'x': 0.0, 'y': 0.0, 'z': 0.0,
                'bond_order_corrected': False, 'skip_reason': 'pdb_block_unparseable',
            })
            continue

        template_smiles = get_smiles(comp_id)
        if not template_smiles:
            out_rows.append({
                'pdb_id': pid, 'ligand_comp_id': comp_id,
                'family': '', 'feature_idx': -1, 'atom_indices_json': '[]',
                'x': 0.0, 'y': 0.0, 'z': 0.0,
                'bond_order_corrected': False, 'skip_reason': 'ccd_smiles_missing',
            })
            continue

        corrected, err = correct_bond_orders(template_smiles, mobile)
        if corrected is None:
            print(f'Stage 4: {pid}:{comp_id} bond-order failed: {err}')
            out_rows.append({
                'pdb_id': pid, 'ligand_comp_id': comp_id,
                'family': '', 'feature_idx': -1, 'atom_indices_json': '[]',
                'x': 0.0, 'y': 0.0, 'z': 0.0,
                'bond_order_corrected': False, 'skip_reason': err,
            })
            continue

        if corrected.GetNumConformers() == 0:
            # The MCS fallback can produce a mol without a conformer. Skip.
            out_rows.append({
                'pdb_id': pid, 'ligand_comp_id': comp_id,
                'family': '', 'feature_idx': -1, 'atom_indices_json': '[]',
                'x': 0.0, 'y': 0.0, 'z': 0.0,
                'bond_order_corrected': True, 'skip_reason': 'no_3d_conformer',
            })
            continue

        conf = corrected.GetConformer()
        feature_idx = 0
        n_emitted = 0
        for fam_code, smarts in smarts_rows:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            for match in corrected.GetSubstructMatches(patt):
                if not match:
                    continue
                xs, ys, zs = [], [], []
                for idx in match:
                    p = conf.GetAtomPosition(idx)
                    xs.append(p.x); ys.append(p.y); zs.append(p.z)
                if not xs:
                    continue
                cx, cy, cz = float(np.mean(xs)), float(np.mean(ys)), float(np.mean(zs))
                nx, ny, nz = apply_4x4_to_point(T, cx, cy, cz)
                out_rows.append({
                    'pdb_id': pid, 'ligand_comp_id': comp_id,
                    'family': FAMILY_NAME[fam_code],
                    'feature_idx': feature_idx,
                    'atom_indices_json': json.dumps(list(match)),
                    'x': nx, 'y': ny, 'z': nz,
                    'bond_order_corrected': True, 'skip_reason': '',
                })
                feature_idx += 1
                n_emitted += 1
        print(f'Stage 4: {pid}:{comp_id} extracted {n_emitted} features')

ligand_features = pd.DataFrame(out_rows)
print(f'Stage 4: emitting {len(ligand_features)} feature rows '
      f'({sum(1 for r in out_rows if r["skip_reason"])} skipped diagnostics)')
