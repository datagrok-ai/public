#name: ProteinLigandInteractionDiagram
#description: Interactive 2D protein-ligand interaction diagram via ProLIF
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, rdkit, mdanalysis, prolif, pdbfixer, openmm, matplotlib, requests]
#input: string protein
#input: string ligand
#input: string ligand_resname
#output: string html

import os
import time
import tempfile
from collections import Counter

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import MDAnalysis as mda
import prolif as plf


# AutoDock atom-type → element mapping (PDBQT col 78-79 → standard PDB col 77-78)
AD_TO_ELEMENT = {
    'A':'C', 'C':'C',
    'N':'N','NA':'N','NS':'N','NX':'N',
    'O':'O','OA':'O','OS':'O',
    'S':'S','SA':'S',
    'H':'H','HD':'H','HS':'H',
    'F':'F','CL':'Cl','BR':'Br','I':'I','P':'P',
    'MG':'Mg','CA':'Ca','MN':'Mn','FE':'Fe','ZN':'Zn','CU':'Cu',
}

SKIP_RESNAMES = {
    'HOH','WAT','H2O','D2O','DOD',
    'NA','CL','K','MG','CA','ZN','MN','FE','CU','NI',
    'SO4','PO4','NO3','ACT','CO3',
    'GOL','EDO','PEG','PG4','DMS','TRS','IMD','BME',
}


def safe_unlink(path, retries=3):
    for _ in range(retries):
        try:
            if os.path.exists(path):
                os.unlink(path)
            return
        except OSError:
            time.sleep(0.1)


def is_pdbqt(text):
    # Ligand PDBQT marker records
    if 'ROOT' in text or 'BRANCH' in text or 'TORSDOF' in text:
        return True
    # Receptor PDBQT: no ROOT/BRANCH but ATOM/HETATM lines carry AutoDock atom
    # types in cols 78-79 (OA, NA, HD, NS, SA, HS, NX, OS) instead of standard
    # element symbols. Sample the first 50 ATOM/HETATM lines.
    ad_specific = {'OA', 'NA', 'HD', 'NS', 'NX', 'OS', 'SA', 'HS'}
    checked = 0
    for line in text.splitlines():
        if line[:6] in ('ATOM  ', 'HETATM') and len(line) >= 79:
            if line[77:79].strip().upper() in ad_specific:
                return True
            checked += 1
            if checked >= 50:
                break
    return False


def pdbqt_to_pdb(text):
    """Strip AutoDock-specific records and restore element symbols in cols 77-78."""
    out = []
    for line in text.splitlines():
        h = line[:6]
        if h.startswith(('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRA', 'TORSDO')):
            continue
        if h in ('ATOM  ', 'HETATM'):
            ad_type = line[77:79].strip().upper() if len(line) >= 78 else ''
            element = AD_TO_ELEMENT.get(ad_type, ad_type[:2].title() if ad_type else '')
            line = line[:66].ljust(76) + element.rjust(2)
        out.append(line)
    return '\n'.join(out)


def shift_negative_residues(pdb_text):
    """Shift residues per chain so all are >= 1."""
    min_per_chain = {}
    for L in pdb_text.splitlines():
        if L[:6] in ('ATOM  ', 'HETATM') and len(L) >= 26:
            try:
                ch, rn = L[21], int(L[22:26])
                if ch not in min_per_chain or rn < min_per_chain[ch]:
                    min_per_chain[ch] = rn
            except ValueError:
                pass
    offsets = {ch: (1 - mn) for ch, mn in min_per_chain.items() if mn < 1}
    if not offsets:
        return pdb_text
    out = []
    for L in pdb_text.splitlines():
        if L[:6] in ('ATOM  ', 'HETATM') and len(L) >= 26:
            try:
                ch, rn = L[21], int(L[22:26])
                if ch in offsets:
                    L = L[:22] + f'{rn + offsets[ch]:>4}' + L[26:]
            except ValueError:
                pass
        out.append(L)
    return '\n'.join(out)


def merge_protein_and_ligand(protein_pdb, ligand_pdb):
    """Concat protein + ligand into one PDB, with TER between."""
    p = [L for L in protein_pdb.splitlines() if not L.startswith('END')]
    l = [L for L in ligand_pdb.splitlines()
         if L[:6] in ('ATOM  ', 'HETATM', 'CONECT')]
    return '\n'.join(p + ['TER'] + l + ['END'])


def detect_ligand_resname(pdb_text):
    counts = Counter(L[17:20].strip() for L in pdb_text.splitlines()
                     if L.startswith('HETATM') and L[17:20].strip() not in SKIP_RESNAMES)
    if not counts:
        raise ValueError('No non-water HETATM ligand found in PDB')
    return counts.most_common(1)[0][0]


def extract_html_from_view(view):
    for getter in ('_repr_html_', 'data', 'to_html'):
        if not hasattr(view, getter):
            continue
        attr = getattr(view, getter)
        try:
            html = attr() if callable(attr) else attr
            if isinstance(html, bytes):
                html = html.decode('utf-8')
            if html:
                return html
        except Exception:
            continue
    raise RuntimeError('Could not extract HTML from LigNetwork view')


# ---------- main ----------

# 1. Normalise protein input (PDBQT → PDB if needed; shift negative residues)
if is_pdbqt(protein):
    protein = pdbqt_to_pdb(protein)
protein = shift_negative_residues(protein)

# 2. Detect ligand resname BEFORE merging — when a separate ligand PDB is
# provided (Docking case), trust IT for the resname rather than the receptor,
# whose HETATMs (buffers, terminal NH3+, ions) would otherwise outvote the
# real ligand via most-common heuristic.
ligand_provided = bool(ligand and ligand.strip())
if not ligand_resname:
    if ligand_provided:
        ligand_resname = detect_ligand_resname(ligand)
    else:
        ligand_resname = detect_ligand_resname(protein)

# 3. Normalise the ligand text (PDBQT → PDB if needed) — keep it SEPARATE for
# now. We deliberately do NOT merge before pdbfixer because pdbfixer drops or
# renames non-standard HETATMs (any ligand whose 3-letter code isn't in its
# template library). Run pdbfixer on the protein alone, then re-attach the
# ligand to the H-fixed PDB afterwards so the ligand survives untouched.
if ligand_provided and is_pdbqt(ligand):
    ligand = pdbqt_to_pdb(ligand)

# 4. pdbfixer for protein hydrogens (on protein-only, no ligand yet)
fd_in, in_path = tempfile.mkstemp(suffix='.pdb'); os.close(fd_in)
fd_out, fixed_path = tempfile.mkstemp(suffix='.pdb'); os.close(fd_out)
fd_combined, combined_path = tempfile.mkstemp(suffix='.pdb'); os.close(fd_combined)
try:
    with open(in_path, 'w') as f:
        f.write(protein)
    fixer = PDBFixer(filename=in_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    with open(fixed_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # 4b. Now merge the ligand into the H-fixed protein. If no separate ligand
    # was provided (BSV case where the cell already contains protein+ligand),
    # the H-fixed PDB already has both — no merge needed.
    if ligand_provided:
        with open(fixed_path) as f:
            fixed_protein = f.read()
        merged = merge_protein_and_ligand(fixed_protein, ligand)
        with open(combined_path, 'w') as f:
            f.write(merged)
        load_path = combined_path
    else:
        load_path = fixed_path

    # 5. MDAnalysis selection (binding-site protein only)
    u = mda.Universe(load_path)
    ligand_ag = u.select_atoms(f'resname {ligand_resname}')
    if len(ligand_ag) == 0:
        raise ValueError(
            f'No atoms with resname {ligand_resname} found in the prepared structure'
        )
    protein_ag = u.select_atoms(
        'protein and byres around 6.0 group ligand', ligand=ligand_ag
    )
    protein_ag.guess_bonds()
    ligand_ag.guess_bonds()

    n_h_lig = len(ligand_ag.select_atoms('element H'))

    # 6. Prepare ligand mol
    if n_h_lig > 0:
        ligand_mol_obj = plf.Molecule.from_mda(ligand_ag)
    else:
        fd_lig, lig_tmp = tempfile.mkstemp(suffix='.pdb'); os.close(fd_lig)
        try:
            ligand_ag.write(lig_tmp)
            lig_rdkit = Chem.MolFromPDBFile(
                lig_tmp, removeHs=False, sanitize=False, proximityBonding=False
            )
            try:
                rdDetermineBonds.DetermineBonds(lig_rdkit, charge=0)
            except Exception:
                lig_rdkit = Chem.MolFromPDBFile(
                    lig_tmp, removeHs=False, sanitize=False, proximityBonding=True
                )
            Chem.SanitizeMol(lig_rdkit)
            # CRITICAL: addResidueInfo=True propagates PDBResidueInfo to new H atoms
            # — without it ProLIF's residue splitting crashes silently on Windows
            lig_rdkit = Chem.AddHs(lig_rdkit, addCoords=True, addResidueInfo=True)
            ligand_mol_obj = plf.Molecule.from_rdkit(lig_rdkit)
        finally:
            safe_unlink(lig_tmp)

    # Protein
    protein_mol_obj = plf.Molecule.from_mda(protein_ag)

    # 7. Fingerprint (n_jobs=1 for Windows safety)
    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand_mol_obj], protein_mol_obj, n_jobs=1)

    # 8. LigNetwork → HTML
    view = fp.plot_lignetwork(ligand_mol_obj, kind='frame', frame=0)
    html = extract_html_from_view(view)
finally:
    safe_unlink(in_path)
    safe_unlink(fixed_path)
    safe_unlink(combined_path)
