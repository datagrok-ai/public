#name: ProteinLigandInteractionDiagram
#description: Interactive 2D protein-ligand interaction diagram via ProLIF
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, rdkit, mdanalysis, prolif, pdbfixer, openmm, matplotlib, requests]
#meta.cache: all
#meta.cache.invalidateOn: 0 0 * * *
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


# ----------------------------------------------------------------------------
# Datagrok pharmacophore SMARTS — used to align ProLIF's interaction detection
# with the Chem package's pharmacophore definitions so a chemist sees the same
# donor/acceptor/charged atoms in both views.
#
# Source of truth: packages/Chem/files/pharmacophore-features.csv
# Combined here per family using recursive SMARTS ($()) so each family becomes
# a single atom expression matching ANY of its constituent patterns.
#
# We only override interactions where SMARTS substitution is well-defined:
#   HBDonor / HBAcceptor / Cationic / Anionic / XBDonor.
# Hydrophobic, PiStacking, CationPi, PiCation, MetalDonor, MetalAcceptor and
# VdWContact keep ProLIF's defaults (geometry-based or different model).
# ----------------------------------------------------------------------------

# Family D — HB donor (donor heavy atom + bonded H). ProLIF's HBDonor expects
# the SMARTS to match a 2-tuple (donor, H), so the trailing [H] is required.
DG_DONOR_SMARTS = (
    '[$([#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]),'
    '$([#8!H0&!$([OH][C,S,P]=O)]),'
    '$([#16!H0]),'
    '$([OX2H1][CX3]=[OX1])]'
    '[H]'
)

# Family A — HB acceptor atom: N/O/F/S acceptors plus aromatic O and nitrile N
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

# Family P — Cationic atom: amines (1°/2°/3°), amidine, guanidine, imidazole,
# any explicit + charge
DG_POSITIVE_SMARTS = (
    '[$([NX3]([CX4])([CX4,#1])[CX4,#1]),'
    '$([CX3](=N)(-N)[!N]),'
    '$(N=[CX3](N)-N),'
    '$([+,+2,+3]),'
    '$([NX3;H2;+0][CX4]),'
    '$([NX3;H1;+0]([CX4])[CX4]),'
    '$(c1c[nH]cn1)]'
)

# Family N — Anionic atom: tetrazole, sulfonate/phosphonate, carboxylate,
# sulfonamide, any explicit − charge
DG_NEGATIVE_SMARTS = (
    '[$(c1nn[nH1]n1),'
    '$([SX4,PX4](=O)(=O)[O-,OH]),'
    '$([CX3,SX3,PX3](=O)[O-,OH]),'
    '$([-,-2,-3]),'
    '$([NR;H1,-1](-C(=O))-[SX4](=O)(=O))]'
)

# Family X — Halogen-bond donor (carbon + bonded halogen). ProLIF's XBDonor
# expects the SMARTS to match a 2-tuple (heavy atom, halogen).
DG_HALOGEN_DONOR_SMARTS = '[#6][Cl,Br,I;X1]'


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


def extract_binding_site_pdb(protein_text, ligand_text, radius=8.0):
    """Filter protein PDB text to only the residues with at least one atom
    within `radius` Å of any ligand atom. Used to shrink the input to
    pdbfixer (which scales linearly with atom count) — typically reduces a
    300-residue protein to 30-50 residues, giving a 5-10× speedup on the
    addMissingHydrogens step.

    Caveat: the cut breaks chain continuity, so the caller must skip
    pdbfixer's findMissingResidues()/addMissingAtoms() — only call
    addMissingHydrogens(). Otherwise pdbfixer would try to "fill" the gaps
    we deliberately removed with phantom amino acids.
    """
    # Collect ligand atom coordinates.
    ligand_xyz = []
    for line in ligand_text.splitlines():
        if line[:6] in ('ATOM  ', 'HETATM') and len(line) >= 54:
            try:
                ligand_xyz.append((
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ))
            except ValueError:
                pass
    if not ligand_xyz:
        return protein_text  # nothing to filter against

    r2 = radius * radius

    # First pass — identify (chain, resid) pairs to keep.
    keep = set()
    for line in protein_text.splitlines():
        if line[:6] not in ('ATOM  ', 'HETATM') or len(line) < 54:
            continue
        try:
            px = float(line[30:38])
            py = float(line[38:46])
            pz = float(line[46:54])
        except ValueError:
            continue
        chain = (line[21] if len(line) > 21 else '').strip() or 'A'
        resid = line[22:26].strip()
        key = (chain, resid)
        if key in keep:
            continue
        for lx, ly, lz in ligand_xyz:
            dx, dy, dz = px - lx, py - ly, pz - lz
            if dx*dx + dy*dy + dz*dz < r2:
                keep.add(key)
                break

    if not keep:
        return protein_text  # safety: keep everything if nothing matched

    # Second pass — output only matching residue lines plus essential headers.
    keep_prefixes = ('HEADER', 'TITLE ', 'CRYST1', 'COMPND', 'MODEL ', 'ENDMDL')
    out = []
    for line in protein_text.splitlines():
        if line[:6] in ('ATOM  ', 'HETATM') and len(line) >= 26:
            chain = (line[21] if len(line) > 21 else '').strip() or 'A'
            resid = line[22:26].strip()
            if (chain, resid) in keep:
                out.append(line)
        elif line[:6] in keep_prefixes or line.startswith('TER') or line.startswith('END'):
            out.append(line)
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

# 2b. Parse the ligand selector (TS may pass "RESNAME CHAIN RESID" to
# disambiguate between multiple instances of the same ligand at different sites)
_parts = ligand_resname.replace(':', ' ').split()
if len(_parts) >= 3:
    sel_resname, sel_chain, sel_resid = _parts[0], _parts[1], _parts[2]
else:
    sel_resname, sel_chain, sel_resid = (_parts[0] if _parts else ligand_resname), None, None

# 3. Normalise the ligand text (PDBQT → PDB if needed) — keep it SEPARATE for
# now. We deliberately do NOT merge before pdbfixer because pdbfixer drops or
# renames non-standard HETATMs (any ligand whose 3-letter code isn't in its
# template library), and OpenMM's PDB writer can also renumber chains/residues.
# Run pdbfixer on the protein alone, then re-attach the ligand to the H-fixed
# PDB afterwards so the ligand survives untouched.
if ligand_provided and is_pdbqt(ligand):
    ligand = pdbqt_to_pdb(ligand)

# 3b. BSV-case fix: if no separate ligand was provided (protein cell contains
# both protein and ligand), extract the matching HETATMs from the protein text
# now and remove them from the protein. This way pdbfixer only sees the
# protein, and we re-merge the original ligand HETATMs (with their original
# chain/resid intact) after pdbfixer.
if not ligand_provided:
    extracted_lines = []
    remaining_lines = []
    for line in protein.splitlines():
        if line.startswith('HETATM') and line[17:20].strip() == sel_resname:
            ln_chain = (line[21] if len(line) > 21 else '').strip() or 'A'
            ln_resid = line[22:26].strip() if len(line) >= 26 else ''
            chain_ok = (sel_chain is None) or (ln_chain == sel_chain)
            resid_ok = (sel_resid is None) or (ln_resid == sel_resid)
            if chain_ok and resid_ok:
                extracted_lines.append(line)
                continue
        remaining_lines.append(line)
    if extracted_lines:
        ligand = '\n'.join(extracted_lines)
        protein = '\n'.join(remaining_lines)
        ligand_provided = True

# 4. pdbfixer for protein hydrogens. Speed optimization: shrink the protein
# to just the binding-site neighborhood (8 Å around the ligand) before running
# pdbfixer. Hydrogenation scales linearly with atom count, so for a 300-residue
# protein this is typically a 5-10× speedup. We only do this when a ligand is
# available to define the neighborhood — for the rare ligand-less case (which
# shouldn't reach this script anyway) the full protein is kept.
fd_in, in_path = tempfile.mkstemp(suffix='.pdb'); os.close(fd_in)
fd_out, fixed_path = tempfile.mkstemp(suffix='.pdb'); os.close(fd_out)
fd_combined, combined_path = tempfile.mkstemp(suffix='.pdb'); os.close(fd_combined)
try:
    if ligand_provided:
        protein_for_fixer = extract_binding_site_pdb(protein, ligand, radius=8.0)
    else:
        protein_for_fixer = protein
    with open(in_path, 'w') as f:
        f.write(protein_for_fixer)

    fixer = PDBFixer(filename=in_path)
    # We only need addMissingHydrogens() — the other pdbfixer steps don't
    # apply to our use case:
    #   findMissingResidues() — detects chain gaps and would try to fill them
    #     with phantom amino acids. After binding-site extraction we have
    #     intentional gaps everywhere; we explicitly do NOT want them filled.
    #     `missingResidues` defaults to {} when this step is skipped.
    #   findMissingAtoms() / addMissingAtoms() — fixes missing side-chain
    #     atoms (real crystal defects). Skipping saves ~100-500ms; in rare
    #     cases an unresolved side-chain atom near the binding site means
    #     slightly approximated H placement, but binding-site residues are
    #     usually well-resolved in deposited structures.
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

    # 5. MDAnalysis selection (binding-site protein only).
    # sel_resname / sel_chain / sel_resid were parsed at step 2b.
    u = mda.Universe(load_path)

    ligand_selection = f'resname {sel_resname}'
    sel_label = sel_resname
    if sel_chain is not None:
        ligand_selection += f' and chainID {sel_chain}'
        sel_label = f'{sel_resname} (chain {sel_chain}, residue {sel_resid})'
    if sel_resid is not None:
        ligand_selection += f' and resid {sel_resid}'

    ligand_ag = u.select_atoms(ligand_selection)
    # Fall back to resname-only if the specific instance can't be found
    # (chain/resid may have been changed by pdbfixer/OpenMM in edge cases).
    if len(ligand_ag) == 0 and (sel_chain is not None or sel_resid is not None):
        ligand_ag = u.select_atoms(f'resname {sel_resname}')
        sel_label = f'{sel_resname} (any instance — original {sel_chain}/{sel_resid} not preserved)'
    if len(ligand_ag) == 0:
        raise ValueError(
            f'No atoms matching {sel_label} found in the prepared structure'
        )
    protein_ag = u.select_atoms(
        # 5.0 Å is enough — ProLIF's longest cutoff is PiStacking at 5.5 Å, but
        # `byres` selects the ENTIRE residue when any of its atoms is within
        # the cutoff, so atoms further out are still included. Tighter than 6.0
        # = smaller atom group = faster MDA→RDKit conversion.
        'protein and byres around 5.0 group ligand', ligand=ligand_ag
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

    # Protein conversion. The MDA→RDKit converter occasionally hits
    # AtomValenceException on residues whose bond perception produces an
    # over-valent atom (modified residues, awkward pdbfixer geometry, glycans,
    # etc.). Try progressively more permissive modes:
    #   1. force=True   — skip strict valence check after bond inference
    #   2. NoImplicit=False — skip bond inference entirely; use implicit Hs.
    #      The molecule's SMARTS matching still works because heavy-atom
    #      connectivity is preserved; only bond orders are simpler.
    try:
        protein_mol_obj = plf.Molecule.from_mda(protein_ag, force=True)
    except Exception as e_force:
        print(f'[ProLIF] from_mda(force=True) failed: {type(e_force).__name__}: {e_force}')
        print('[ProLIF] Retrying with NoImplicit=False (skip bond-order inference)')
        try:
            protein_mol_obj = plf.Molecule.from_mda(protein_ag, NoImplicit=False)
        except TypeError:
            # Very old MDA without these kwargs — last-ditch plain call
            protein_mol_obj = plf.Molecule.from_mda(protein_ag)

    # 7. Fingerprint with Datagrok-aligned SMARTS where it makes sense.
    # If our combined SMARTS fail to parse for any reason, fall back to
    # ProLIF's defaults so the script still produces a useful result.
    interactions = [
        'Hydrophobic', 'HBDonor', 'HBAcceptor',
        'PiStacking', 'Cationic', 'Anionic', 'CationPi', 'PiCation',
        'XBDonor', 'XBAcceptor', 'MetalDonor', 'MetalAcceptor', 'VdWContact',
    ]
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
        # Datagrok SMARTS rejected by ProLIF/RDKit — fall back so we still get a diagram
        print(f'[ProLIF] Datagrok-aligned SMARTS rejected ({e}); using ProLIF defaults')
        fp = plf.Fingerprint(interactions=interactions)
    fp.run_from_iterable([ligand_mol_obj], protein_mol_obj, n_jobs=1)

    # 8. LigNetwork → HTML
    view = fp.plot_lignetwork(ligand_mol_obj, kind='frame', frame=0)
    html = extract_html_from_view(view)
finally:
    safe_unlink(in_path)
    safe_unlink(fixed_path)
    safe_unlink(combined_path)
