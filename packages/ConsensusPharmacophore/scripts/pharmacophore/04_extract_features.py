#name: PharmacophoreStage4ExtractFeatures
#description: Stage 4 - per-ligand pharmacophore extraction via ProLIF. Each detected protein-ligand interaction emits one feature row; ligand-side atom centroid lifted into the consensus frame.
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, rdkit, mdanalysis, prolif, pdbfixer, openmm, matplotlib, requests, pip, {pip: [openbabel-wheel]}]
#meta.cache: none
#input: dataframe aligned_structures {caption: Stage 2b pass-2 aligned structures}
#input: dataframe pocket_atoms {caption: Stage 3 pocket atoms (used to discover ligand comp_ids per PDB)}
#input: string selected_pdb_ids = "[]" {caption: JSON array of PDB IDs to keep}
#output: dataframe ligand_features {caption: One row per ProLIF interaction}

# ProLIF-based pharmacophore extractor for the Consensus Pharmacophore pipeline.
#
# Replaces the earlier SMARTS-only Stage 4 (archived at
# files/_reference/04_extract_features_smarts_only.py). The earlier version
# matched the Datagrok pharmacophore SMARTS dictionary against the ligand atoms
# in isolation: it had no awareness of which features actually contact the
# protein, so a ligand N-H pointing into bulk solvent contributed the same
# Donor feature as one H-bonded to a hinge backbone carbonyl.
#
# This version runs ProLIF over the (protein + ligand) complex and emits a
# row only when ProLIF detects a productive interaction. The ligand atom
# centroid of each detected interaction is the consensus pharmacophore
# point candidate; Stage 5a then clusters those candidates per family.
#
# Implementation notes:
#   - Environment string is COPIED VERBATIM from BiostructureViewer's ProLIF
#     script so the conda env (which is MD5-keyed by the env line) is shared
#     across both packages. Saves a full env build on the worker.
#   - All the PDB parsing helpers (PDBQT detection, residue shifting, binding-
#     site extraction, protein/ligand merging, atom-name uniquification) are
#     copied / adapted from `BiostructureViewer/scripts/protein_ligand_interactions.py`
#     because Datagrok scripts can't share helper modules across packages.
#   - We loop over rows of `aligned_structures` internally, one ProLIF run per
#     PDB. ProLIF is single-threaded; sequential loop matches BSV's design.
#   - The ligand atom indices on a ProLIF interaction reference the
#     plf.Molecule (post-pdbfixer), which is renumbered relative to the
#     original PDB. We use plf.Molecule positions directly — pdbfixer only
#     ADDS Hs/missing atoms; existing heavy atoms keep their original
#     coordinates, which is what `transform_4x4_json` was computed against.

import json
import logging
import os
import tempfile
from collections import Counter
from io import StringIO
from itertools import islice

# certifi → SSL_CERT_FILE so pdbfixer's CCD template fetch from RCSB works
# in containers without a system CA bundle. Matches the BSV script's
# rationale verbatim: without this, `addMissingHydrogens` silently fails
# to add Hs to non-standard ligands and downstream bond perception breaks.
try:
    import certifi
    os.environ.setdefault('SSL_CERT_FILE', certifi.where())
    os.environ.setdefault('REQUESTS_CA_BUNDLE', certifi.where())
except ImportError:
    pass

import numpy as np
import pandas as pd
from rdkit import Chem
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import MDAnalysis as mda
from MDAnalysis.lib.util import NamedStream
from MDAnalysis.guesser.tables import vdwradii as MDA_VDWRADII
import prolif as plf


logger = logging.getLogger('pharmacophore.stage4.prolif')


# Mixed-case VdW radii aliases — see BSV's protein_ligand_interactions.py
# for the full story. Required so guess_bonds() on ligands with halogens
# doesn't raise "vdw radii for types: Cl".
_VDW_RADII = dict(MDA_VDWRADII)
for _src, _aliases in (
    ('CL', ('Cl', 'CL')),
    ('BR', ('Br', 'BR')),
    ('I',  ('I', 'i', 'IOD')),
    ('NA', ('Na', 'NA')),
    ('MG', ('Mg', 'MG')),
    ('K',  ('K', 'k')),
    ('CA', ('Ca', 'CA')),
    ('ZN', ('Zn', 'ZN')),
    ('FE', ('Fe', 'FE')),
):
    if _src in _VDW_RADII:
        for _a in _aliases:
            _VDW_RADII.setdefault(_a, _VDW_RADII[_src])


# ----------------------------------------------------------------------------
# Datagrok pharmacophore SMARTS — copied verbatim from BSV's
# protein_ligand_interactions.py so ProLIF flags the same donor/acceptor /
# charged atoms the rest of the platform shows.
# ----------------------------------------------------------------------------

DG_DONOR_SMARTS = (
    '[$([#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]),'
    '$([#8!H0&!$([OH][C,S,P]=O)]),'
    '$([#16!H0]),'
    '$([#7;+;!H0]),'
    '$([OX2H1][CX3]=[OX1])]'
    '[H]'
)
DG_ACCEPTOR_SMARTS = (
    '[$([#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)&!$([#7;+])]),'
    '$([O;!$([OX2](C)C=O);!$(*(~a)~a)]),'
    '$([F;$(F-[#6])]),'
    '$([SX2;H0;v2]),'
    '$([SX2;H1;+0]),'
    '$([S;-]),'
    '$([o+0]),'
    '$([NX1]#[CX2])]'
)
DG_POSITIVE_SMARTS = (
    '[$([NX3]([CX4])([CX4,#1])[CX4,#1])&!$([NX3]-*=[!#6]),'
    '$([CX3](=N)(-N)[!N]),'
    '$(N=[CX3](N)-N),'
    '$([+,+2,+3])&!$(*[-,-2,-3]),'
    '$([NX3;H2;+0][CX4])&!$([NX3;H2](C)C=[O,N,S]),'
    '$([NX3;H1;+0]([CX4])[CX4])&!$([NX3;H1](C)(C)C=[O,N,S]),'
    '$(c1c[nH]cn1)]'
)
DG_NEGATIVE_SMARTS = (
    '[$(c1nn[nH1]n1),'
    '$([SX4,PX4](=O)(=O)[O-,OH]),'
    '$([CX3,SX3,PX3](=O)[O-,OH]),'
    '$([-,-2,-3])&!$(*[+,+2,+3]),'
    '$([NR;H1,-1](-C(=O))-[SX4](=O)(=O))]'
)
DG_HALOGEN_DONOR_SMARTS = '[#6][Cl,Br,I;X1]'


# ----------------------------------------------------------------------------
# Interaction → Datagrok family. ProLIF detects more interaction types than
# the 7 families we render; we map the ones that fit and skip the rest.
#
# `family_name` is the human-readable name in family-map.ts (Donor/Acceptor/…).
# Skipped: VdWContact (geometric fallback, way too many hits — would mostly
#          drown out specific contacts), MetalDonor/MetalAcceptor (rare, no
#          family in the v1 taxonomy).
# ----------------------------------------------------------------------------

INT_TO_FAMILY = {
    'HBDonor':    ('D', 'Donor'),
    'HBAcceptor': ('A', 'Acceptor'),
    'PiStacking': ('a', 'Aromatic'),   # ligand has aromatic ring
    'PiCation':   ('a', 'Aromatic'),   # ligand aromatic, protein cation
    'Hydrophobic': ('H', 'Hydrophobic'),
    'Cationic':   ('P', 'Positive'),
    'CationPi':   ('P', 'Positive'),   # ligand cation, protein aromatic
    'Anionic':    ('N', 'Negative'),
    'XBDonor':    ('X', 'Halogen'),
}
SKIPPED_INTERACTIONS = frozenset({
    'VdWContact', 'MetalDonor', 'MetalAcceptor', 'XBAcceptor',
})


# AutoDock atom-type → element mapping. Carried over from BSV for the rare
# case where Stage 1 received a PDBQT instead of a PDB.
AD_TO_ELEMENT = {
    'A': 'C', 'C': 'C',
    'N': 'N', 'NA': 'N', 'NS': 'N', 'NX': 'N',
    'O': 'O', 'OA': 'O', 'OS': 'O',
    'S': 'S', 'SA': 'S',
    'H': 'H', 'HD': 'H', 'HS': 'H',
    'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I', 'P': 'P',
    'MG': 'Mg', 'CA': 'Ca', 'MN': 'Mn', 'FE': 'Fe', 'ZN': 'Zn', 'CU': 'Cu',
}

ATOM_HETATM_PREFIXES = ('ATOM  ', 'HETATM')
PDB_HEADER_PREFIXES = ('HEADER', 'TITLE ', 'CRYST1', 'COMPND', 'MODEL ', 'ENDMDL')
AUTODOCK_ATOM_TYPES = frozenset({'OA', 'NA', 'HD', 'NS', 'NX', 'OS', 'SA', 'HS'})
PDBQT_LIGAND_MARKERS = ('ROOT', 'BRANCH', 'TORSDOF')
PDBQT_RECORDS_TO_DROP = ('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRA', 'TORSDO')
MERGE_LIGAND_RECORDS = ('ATOM  ', 'HETATM', 'CONECT')

# Defensive filter: standard amino acids should never appear in `ligand_by_pdb`.
# Stage 3 is supposed to mark only HETATM-record drug ligands as `ligand_seed=True`,
# but the Datagrok DataFrame round-trip has been observed coercing a few pocket-Cα
# rows' booleans into truthy values — feeding ALA/ARG/ASP/CYS into our ProLIF
# loop and emitting useless `no HETATM with resname X in PDB` diagnostic rows.
# Skipping these resnames here is harmless even when Stage 3 is correct: a real
# small-molecule drug never has a standard-AA 3-letter code.
STANDARD_AA_RESNAMES = frozenset({
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'SEC', 'PYL',
    'HOH', 'WAT', 'H2O', 'D2O', 'DOD',
})


def _is_truthy_ligand_seed(v):
    """Robust truthy check for the `ligand_seed` cell value.
    Handles Python bool, numpy bool, 0/1, and the 'False'/'False' string form
    that the Datagrok DataFrame round-trip can produce (where `bool('False')`
    would otherwise return True because any non-empty string is truthy)."""
    if v is True or v == 1:
        return True
    if isinstance(v, str):
        return v.strip().lower() in ('true', '1', 'yes')
    return False


def is_pdbqt_lines(lines):
    for line in lines:
        if line.startswith(PDBQT_LIGAND_MARKERS):
            return True
        if line.startswith(ATOM_HETATM_PREFIXES):
            break
    candidates = (line for line in lines
                  if line.startswith(ATOM_HETATM_PREFIXES) and len(line) >= 79)
    return any(line[77:79].strip().upper() in AUTODOCK_ATOM_TYPES
               for line in islice(candidates, 50))


def pdbqt_to_pdb_lines(lines):
    out = []
    for line in lines:
        if line.startswith(PDBQT_RECORDS_TO_DROP):
            continue
        if line.startswith(ATOM_HETATM_PREFIXES):
            ad_type = line[77:79].strip().upper() if len(line) >= 79 else ''
            element = AD_TO_ELEMENT.get(ad_type, ad_type[:2].title() if ad_type else '')
            line = line[:66].ljust(76) + element.rjust(2)
        out.append(line)
    return out


def shift_negative_residues_lines(lines):
    """Shift residues per chain so all are >= 1 (matches BSV)."""
    parsed = []
    min_per_chain = {}
    for idx, L in enumerate(lines):
        if L.startswith(ATOM_HETATM_PREFIXES) and len(L) >= 26:
            try:
                ch, rn = L[21], int(L[22:26])
            except ValueError:
                continue
            parsed.append((idx, ch, rn))
            if ch not in min_per_chain or rn < min_per_chain[ch]:
                min_per_chain[ch] = rn
    offsets = {ch: (1 - mn) for ch, mn in min_per_chain.items() if mn < 1}
    if not offsets:
        return lines
    out = list(lines)
    for idx, ch, rn in parsed:
        if ch in offsets:
            L = out[idx]
            out[idx] = L[:22] + f'{rn + offsets[ch]:>4}' + L[26:]
    return out


def extract_ligand_lines_for_resname(lines, comp_id):
    """Return the HETATM lines whose resname column equals comp_id (upper)."""
    cid = comp_id.upper()
    out = []
    for line in lines:
        if line.startswith('HETATM') and len(line) >= 20 \
                and line[17:20].strip().upper() == cid:
            out.append(line)
    return out


def remove_lines_for_resname(lines, comp_id):
    """Return `lines` with HETATM entries for `comp_id` removed."""
    cid = comp_id.upper()
    out = []
    for line in lines:
        if line.startswith('HETATM') and len(line) >= 20 \
                and line[17:20].strip().upper() == cid:
            continue
        out.append(line)
    return out


def extract_binding_site_lines(protein_lines, ligand_lines, radius=8.0):
    """Filter protein lines to residues within radius of any ligand atom.
    Speed-only optimization for pdbfixer (it scales linearly with atom count).
    """
    ligand_xyz = []
    for line in ligand_lines:
        if line.startswith(ATOM_HETATM_PREFIXES) and len(line) >= 54:
            try:
                ligand_xyz.append((
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ))
            except ValueError:
                pass
    if not ligand_xyz:
        return protein_lines
    r2 = radius * radius
    chains, resids, xs, ys, zs = [], [], [], [], []
    for line in protein_lines:
        if not line.startswith(ATOM_HETATM_PREFIXES) or len(line) < 54:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        chains.append((line[21] if len(line) > 21 else '').strip() or 'A')
        resids.append(line[22:26].strip())
        xs.append(x); ys.append(y); zs.append(z)
    if not xs:
        return protein_lines
    prot_xyz = np.column_stack((xs, ys, zs))
    lig_xyz_np = np.asarray(ligand_xyz, dtype=float)
    diff = prot_xyz[:, None, :] - lig_xyz_np[None, :, :]
    min_d2 = (diff * diff).sum(axis=-1).min(axis=1)
    in_range = np.flatnonzero(min_d2 < r2)
    keep = {(chains[i], resids[i]) for i in in_range}
    if not keep:
        return protein_lines
    out = []
    for line in protein_lines:
        if line.startswith(ATOM_HETATM_PREFIXES) and len(line) >= 26:
            chain = (line[21] if len(line) > 21 else '').strip() or 'A'
            resid = line[22:26].strip()
            if (chain, resid) in keep:
                out.append(line)
        elif line.startswith(PDB_HEADER_PREFIXES) or line.startswith(('TER', 'END')):
            out.append(line)
    return out


def merge_protein_and_ligand_lines(protein_lines, ligand_lines):
    """Concat protein + ligand with TER between. Renumber ligand serials so
    they don't collide with the protein's, and patch CONECT references.
    """
    max_serial = 0
    p_kept = []
    for L in protein_lines:
        if L.startswith('END'):
            continue
        p_kept.append(L)
        if L.startswith(ATOM_HETATM_PREFIXES) and len(L) >= 11:
            try:
                max_serial = max(max_serial, int(L[6:11]))
            except ValueError:
                pass
    raw_ligand = [L for L in ligand_lines if L.startswith(MERGE_LIGAND_RECORDS)]
    remap = {}
    next_serial = max_serial + 1
    for L in raw_ligand:
        if L.startswith(ATOM_HETATM_PREFIXES) and len(L) >= 11:
            try:
                old = int(L[6:11])
            except ValueError:
                continue
            if old not in remap:
                remap[old] = next_serial
                next_serial += 1
    l_kept = []
    for L in raw_ligand:
        if L.startswith(ATOM_HETATM_PREFIXES) and len(L) >= 11:
            try:
                old = int(L[6:11])
                if old in remap:
                    L = L[:6] + f'{remap[old]:>5}' + L[11:]
            except ValueError:
                pass
            l_kept.append(L)
        elif L.startswith('CONECT'):
            try:
                serials = []
                pos = 6
                while pos + 5 <= len(L):
                    chunk = L[pos:pos + 5].strip()
                    if not chunk:
                        break
                    serials.append(int(chunk))
                    pos += 5
                if serials:
                    new_serials = [remap.get(s, s) for s in serials]
                    L = 'CONECT' + ''.join(f'{s:>5}' for s in new_serials)
            except ValueError:
                continue
            l_kept.append(L)
    return p_kept + ['TER'] + l_kept + ['END']


def ensure_unique_atom_names(lines):
    """Rewrite atom names within each (chain, resid) so every atom has a
    unique name. AutoDock poses use generic 'C', 'N', 'O' names that collide
    under pdbfixer's (residue, atom_name) dedup. Idempotent on unique input.
    """
    counters = {}
    out = []
    for L in lines:
        if not L.startswith(ATOM_HETATM_PREFIXES) or len(L) < 27:
            out.append(L)
            continue
        chain = L[21]
        resid = L[22:26]
        elem = L[76:78].strip() if len(L) >= 78 else ''
        if not elem:
            elem = L[12:16].strip().rstrip('0123456789') or 'X'
        key = (chain, resid)
        n = counters.get(key, {}).get(elem, 0) + 1
        counters.setdefault(key, {})[elem] = n
        new_name = f'{elem}{n}'.ljust(4)[:4]
        out.append(L[:12] + new_name + L[16:])
    return out


def apply_4x4_to_point(T, x, y, z):
    """T is a 4×4 numpy array (row-major from Stage 2b); returns (nx, ny, nz)."""
    R = T[:3, :3]
    t = T[:3, 3]
    p = np.array([x, y, z], dtype=float)
    out = R @ p + t
    return float(out[0]), float(out[1]), float(out[2])


# ----------------------------------------------------------------------------
# Per-PDB ProLIF processing
# ----------------------------------------------------------------------------

def _build_ligand_mol_obj(ligand_ag):
    """Hs-present → from_mda (fast). Hs-absent → OpenBabel rebuild → from_rdkit
    (slow but robust). Mirrors BSV's two-tier setup."""
    n_h_lig = len(ligand_ag.select_atoms('element H'))
    if n_h_lig > 0:
        try:
            return plf.Molecule.from_mda(ligand_ag)
        except Exception as e_mda:
            logger.warning('from_mda failed (%s: %s); falling back to OpenBabel',
                           type(e_mda).__name__, e_mda)
            ligand_ag = ligand_ag.select_atoms('not element H')

    # OpenBabel rebuild path
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as _tmp:
        ligpath = _tmp.name
    lig_pdb_text = ''
    try:
        ligand_ag.write(ligpath)
        with open(ligpath, 'r') as f:
            lig_pdb_text = f.read()
        from openbabel import pybel as _pybel
        ob_mol = next(_pybel.readfile('pdb', ligpath))
        ob_mol.OBMol.PerceiveBondOrders()
        ob_mol.addh()
        mol_block = ob_mol.write('mol')
    finally:
        try:
            os.unlink(ligpath)
        except OSError:
            pass

    lig_rdkit = None
    try:
        lig_rdkit = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
    except Exception:
        pass
    if lig_rdkit is None:
        try:
            m = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
            if m is not None:
                Chem.SanitizeMol(m, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
                lig_rdkit = m
        except Exception:
            pass
    if lig_rdkit is None:
        try:
            m = Chem.MolFromPDBBlock(lig_pdb_text, removeHs=False, sanitize=False,
                                     proximityBonding=True)
            if m is not None:
                try:
                    Chem.SanitizeMol(m)
                except Exception:
                    Chem.SanitizeMol(m, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
                lig_rdkit = Chem.AddHs(m, addCoords=True, addResidueInfo=True)
        except Exception:
            pass
    if lig_rdkit is None:
        raise RuntimeError('All ligand parse paths failed (MolBlock strict/lenient, '
                           'PDB proximity-bonding).')
    return plf.Molecule.from_rdkit(lig_rdkit)


def _ligand_atom_positions(ligand_mol_obj):
    """Return Nx3 numpy array of ALL ligand atom positions in the same frame
    as the original PDB (pdbfixer doesn't move existing heavy atoms)."""
    rdkit_mol = ligand_mol_obj  # plf.Molecule inherits rdkit.Chem.Mol
    conf = rdkit_mol.GetConformer()
    n = rdkit_mol.GetNumAtoms()
    pos = np.empty((n, 3), dtype=float)
    for i in range(n):
        p = conf.GetAtomPosition(i)
        pos[i, 0] = p.x; pos[i, 1] = p.y; pos[i, 2] = p.z
    return pos


def _run_prolif_for_pdb(pdb_text, comp_id, T):
    """Run ProLIF on a single (PDB text, ligand comp_id) pair and return a
    list of dict rows (one per detected interaction, ligand-centroid lifted
    into consensus frame via T).
    """
    protein_lines = pdb_text.splitlines()
    if is_pdbqt_lines(protein_lines):
        protein_lines = pdbqt_to_pdb_lines(protein_lines)
    protein_lines = shift_negative_residues_lines(protein_lines)

    ligand_lines = extract_ligand_lines_for_resname(protein_lines, comp_id)
    if not ligand_lines:
        return [{'skip_reason': f'no HETATM with resname {comp_id} in PDB',
                 'family': '', 'interaction_type': '', 'residue': '',
                 'x': 0.0, 'y': 0.0, 'z': 0.0,
                 'distance': 0.0, 'atom_indices_json': '[]'}]

    protein_only_lines = remove_lines_for_resname(protein_lines, comp_id)
    ligand_lines = ensure_unique_atom_names(ligand_lines)

    # Shrink the protein to the binding-site neighbourhood (8 Å) — 5-10×
    # speedup on pdbfixer for a typical 300-residue protein.
    protein_for_fixer_lines = extract_binding_site_lines(
        protein_only_lines, ligand_lines, radius=8.0)
    structure_for_fixer_lines = merge_protein_and_ligand_lines(
        protein_for_fixer_lines, ligand_lines)
    structure_for_fixer = '\n'.join(structure_for_fixer_lines)

    fixer = PDBFixer(pdbfile=StringIO(structure_for_fixer))
    fixer.missingResidues = {}
    try:
        fixer.findNonstandardResidues()
        if fixer.nonstandardResidues:
            fixer.replaceNonstandardResidues()
    except Exception as e:
        logger.warning('pdbfixer replaceNonstandardResidues failed: %s', e)
    try:
        fixer.findMissingAtoms()
        if fixer.missingAtoms or fixer.missingTerminals:
            fixer.addMissingAtoms()
    except Exception as e:
        logger.warning('pdbfixer addMissingAtoms failed: %s', e)
    fixer.addMissingHydrogens(pH=7.0)

    buf = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, buf)
    merged_pdb = buf.getvalue()

    u = mda.Universe(NamedStream(StringIO(merged_pdb), 'merged.pdb'))
    ligand_ag = u.select_atoms(f'resname {comp_id}')
    if len(ligand_ag) == 0:
        return [{'skip_reason': f'ligand {comp_id} disappeared during pdbfixer',
                 'family': '', 'interaction_type': '', 'residue': '',
                 'x': 0.0, 'y': 0.0, 'z': 0.0,
                 'distance': 0.0, 'atom_indices_json': '[]'}]

    # 6.5 Å: covers ProLIF's longest interaction cutoff (EdgeToFace pi-stacking).
    protein_ag = u.select_atoms(
        'protein and byres around 6.5 group ligand', ligand=ligand_ag)
    protein_ag.guess_bonds(vdwradii=_VDW_RADII)
    ligand_ag.guess_bonds(vdwradii=_VDW_RADII)

    ligand_mol_obj = _build_ligand_mol_obj(ligand_ag)
    protein_mol_obj = plf.Molecule.from_mda(protein_ag)

    interactions = list(INT_TO_FAMILY.keys())  # only the families we render
    parameters = {
        'HBDonor':    {'donor': DG_DONOR_SMARTS, 'acceptor': DG_ACCEPTOR_SMARTS},
        'HBAcceptor': {'donor': DG_DONOR_SMARTS, 'acceptor': DG_ACCEPTOR_SMARTS},
        'Cationic':   {'cation': DG_POSITIVE_SMARTS, 'anion': DG_NEGATIVE_SMARTS},
        'Anionic':    {'cation': DG_POSITIVE_SMARTS, 'anion': DG_NEGATIVE_SMARTS},
        'XBDonor':    {'donor': DG_HALOGEN_DONOR_SMARTS},
    }
    try:
        fp = plf.Fingerprint(interactions=interactions, parameters=parameters)
    except Exception as e:
        logger.warning('Datagrok-aligned SMARTS rejected (%s); using ProLIF defaults', e)
        fp = plf.Fingerprint(interactions=interactions)
    fp.run_from_iterable([ligand_mol_obj], protein_mol_obj, n_jobs=1)

    lig_positions = _ligand_atom_positions(ligand_mol_obj)

    rows = []
    # ifp[frame_id] = { (lig_resid, prot_resid): {interaction: tuple[metadata, ...]} }
    # We only ran one frame (frame_id = 0).
    try:
        frame_ifp = fp.ifp[0]
    except (KeyError, IndexError, TypeError):
        return [{'skip_reason': 'fp.ifp[0] missing — ProLIF produced no frames',
                 'family': '', 'interaction_type': '', 'residue': '',
                 'x': 0.0, 'y': 0.0, 'z': 0.0,
                 'distance': 0.0, 'atom_indices_json': '[]'}]

    for (lig_resid, prot_resid), interactions_dict in frame_ifp.items():
        residue_str = str(prot_resid)  # e.g. 'TYR123.A'
        for int_name, metadata_tuple in interactions_dict.items():
            if int_name in SKIPPED_INTERACTIONS:
                continue
            family_code_name = INT_TO_FAMILY.get(int_name)
            if family_code_name is None:
                continue
            _fam_code, fam_name = family_code_name
            for md in metadata_tuple:
                # ProLIF metadata keys vary by version: try the typed
                # `parent_indices['ligand']` first (post-v1.0), fall back to
                # `indices['ligand']` for older releases. Both reference the
                # ligand mol's atom numbering.
                lig_idx = None
                indices = md.get('parent_indices') or md.get('indices') or {}
                lig_idx = indices.get('ligand')
                if lig_idx is None or len(lig_idx) == 0:
                    continue
                # Centroid the ligand-side atoms; lift into consensus frame.
                xs = lig_positions[list(lig_idx), 0]
                ys = lig_positions[list(lig_idx), 1]
                zs = lig_positions[list(lig_idx), 2]
                cx = float(np.mean(xs)); cy = float(np.mean(ys)); cz = float(np.mean(zs))
                nx, ny, nz = apply_4x4_to_point(T, cx, cy, cz)
                rows.append({
                    'family': fam_name,
                    'interaction_type': int_name,
                    'residue': residue_str,
                    'distance': float(md.get('distance', 0.0) or 0.0),
                    'atom_indices_json': json.dumps(list(int(i) for i in lig_idx)),
                    'x': nx, 'y': ny, 'z': nz,
                    'skip_reason': '',
                })
    return rows


# ----------------------------------------------------------------------------
# Main body
# ----------------------------------------------------------------------------

if not isinstance(aligned_structures, pd.DataFrame) or len(aligned_structures) == 0:
    raise ValueError('Stage 4: aligned_structures is empty.')
if not isinstance(pocket_atoms, pd.DataFrame) or len(pocket_atoms) == 0:
    raise ValueError('Stage 4: pocket_atoms is empty.')

try:
    selected = set(json.loads(selected_pdb_ids or '[]'))
except (json.JSONDecodeError, TypeError):
    selected = set()
use_all = (not selected) or ('all' in selected)

# Per-PDB ligand comp_ids from Stage 3's `ligand_seed=True` flag.
ligand_by_pdb = {}
skipped_aa_count = 0
for i in range(len(pocket_atoms)):
    row = pocket_atoms.iloc[i]
    if not _is_truthy_ligand_seed(row.get('ligand_seed')):
        continue
    pid = str(row['pdb_id'])
    comp = str(row.get('res_name', '')).strip().upper()
    if not comp:
        continue
    if comp in STANDARD_AA_RESNAMES:
        skipped_aa_count += 1
        continue
    ligand_by_pdb.setdefault(pid, set()).add(comp)
if skipped_aa_count:
    print(f'Stage 4 ProLIF: skipped {skipped_aa_count} pocket_atoms rows whose '
          f'res_name was a standard amino acid (defensive filter — Stage 3 should '
          f'have marked these as ligand_seed=False; see STANDARD_AA_RESNAMES).')
print(f'Stage 4 ProLIF: ligands to process: '
      f'{ {p: sorted(c) for p, c in ligand_by_pdb.items()} }')

out_rows = []
n_pdb_total = 0
n_interactions_total = 0
for i in range(len(aligned_structures)):
    src = aligned_structures.iloc[i].to_dict()
    pid = str(src['pdb_id'])
    if not use_all and pid not in selected:
        continue
    if pid not in ligand_by_pdb:
        print(f'Stage 4 ProLIF: {pid} has no ligand in pocket_atoms — skip')
        continue

    pdb_text = str(src.get('original_pdb', '') or '')
    if not pdb_text:
        out_rows.append({
            'pdb_id': pid, 'ligand_comp_id': '',
            'family': '', 'interaction_type': '', 'residue': '',
            'distance': 0.0, 'atom_indices_json': '[]',
            'x': 0.0, 'y': 0.0, 'z': 0.0,
            'skip_reason': 'pdb_text_missing',
        })
        continue

    try:
        T = np.array(json.loads(src['transform_4x4_json']), dtype=float).reshape(4, 4)
    except (ValueError, KeyError, json.JSONDecodeError):
        T = np.eye(4)

    n_pdb_total += 1
    for comp_id in sorted(ligand_by_pdb[pid]):
        try:
            rows = _run_prolif_for_pdb(pdb_text, comp_id, T)
        except Exception as e:
            print(f'Stage 4 ProLIF: {pid}:{comp_id} failed — {type(e).__name__}: {e}')
            rows = [{'family': '', 'interaction_type': '', 'residue': '',
                     'distance': 0.0, 'atom_indices_json': '[]',
                     'x': 0.0, 'y': 0.0, 'z': 0.0,
                     'skip_reason': f'{type(e).__name__}: {e}'}]
        n_good = sum(1 for r in rows if not r.get('skip_reason'))
        n_interactions_total += n_good
        for r in rows:
            r['pdb_id'] = pid
            r['ligand_comp_id'] = comp_id
            out_rows.append(r)
        print(f'Stage 4 ProLIF: {pid}:{comp_id} → {n_good} interactions')

if out_rows:
    ligand_features = pd.DataFrame(out_rows)
    # Column ordering — pdb_id / ligand_comp_id first, then family/interaction/residue,
    # then coords. Keeps the table readable when previewed in Datagrok's grid.
    preferred = ['pdb_id', 'ligand_comp_id', 'family', 'interaction_type', 'residue',
                 'distance', 'atom_indices_json', 'x', 'y', 'z', 'skip_reason']
    keep = [c for c in preferred if c in ligand_features.columns]
    ligand_features = ligand_features[keep]
else:
    ligand_features = pd.DataFrame(columns=[
        'pdb_id', 'ligand_comp_id', 'family', 'interaction_type', 'residue',
        'distance', 'atom_indices_json', 'x', 'y', 'z', 'skip_reason',
    ])

print(f'Stage 4 ProLIF: emitting {len(ligand_features)} rows '
      f'({n_interactions_total} real interactions across {n_pdb_total} PDBs; '
      f'{sum(1 for r in out_rows if r.get("skip_reason"))} skipped diagnostics)')
