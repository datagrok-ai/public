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

import logging
from collections import Counter
from io import StringIO

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import MDAnalysis as mda
from MDAnalysis.lib.util import NamedStream
import prolif as plf

logger = logging.getLogger('prolif')


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

SKIP_RESNAMES = frozenset({
    'HOH','WAT','H2O','D2O','DOD',
    'NA','CL','K','MG','CA','ZN','MN','FE','CU','NI',
    'SO4','PO4','NO3','ACT','CO3',
    'GOL','EDO','PEG','PG4','DMS','TRS','IMD','BME',
})

# Hoisted constants — defined at module scope so they're built once at import
# rather than recreated on every helper call (and, for some, every loop iter).

# PDB record prefixes used in fixed-column line parsing
ATOM_HETATM_PREFIXES = ('ATOM  ', 'HETATM')
PDB_HEADER_PREFIXES = ('HEADER', 'TITLE ', 'CRYST1', 'COMPND', 'MODEL ', 'ENDMDL')

# AutoDock receptor PDBQTs lack ROOT/BRANCH but carry these atom-type codes
# in cols 78-79; sample-based detection looks for them.
AUTODOCK_ATOM_TYPES = frozenset({'OA', 'NA', 'HD', 'NS', 'NX', 'OS', 'SA', 'HS'})

# AutoDock ligand PDBQTs carry these record types alongside ATOM/HETATM
PDBQT_LIGAND_MARKERS = ('ROOT', 'BRANCH', 'TORSDOF')
# All AutoDock-specific record prefixes that should be dropped during PDBQT→PDB
PDBQT_RECORDS_TO_DROP = ('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRA', 'TORSDO')

# Records to preserve when concatenating ligand HETATMs into a fixed protein
MERGE_LIGAND_RECORDS = ('ATOM  ', 'HETATM', 'CONECT')


# Helper functions take/return list[str] (one element per PDB line) so the
# main flow only calls splitlines() once per text input. Calling splitlines on
# a few-thousand-line PDB is ~1-5ms each — small individually, but it adds up
# across 6-8 calls per row. Operating on lists everywhere makes the cost
# proportional to the work, not to the number of helpers.

def is_pdbqt_lines(lines):
    # Ligand PDBQT marker records — return early if found
    for line in lines:
        if line.startswith(PDBQT_LIGAND_MARKERS):
            return True
        if line.startswith(ATOM_HETATM_PREFIXES):
            break
    # Receptor PDBQT: no ROOT/BRANCH, but ATOM/HETATM lines carry AutoDock atom
    # types in cols 78-79 (OA, NA, HD, NS, SA, HS, NX, OS) instead of standard
    # element symbols. Sample the first 50 ATOM/HETATM lines.
    checked = 0
    for line in lines:
        if line[:6] in ATOM_HETATM_PREFIXES and len(line) >= 79:
            if line[77:79].strip().upper() in AUTODOCK_ATOM_TYPES:
                return True
            checked += 1
            if checked >= 50:
                break
    return False


def pdbqt_to_pdb_lines(lines):
    """Strip AutoDock-specific records and restore element symbols in cols 77-78."""
    out = []
    for line in lines:
        h = line[:6]
        if h.startswith(PDBQT_RECORDS_TO_DROP):
            continue
        if h in ATOM_HETATM_PREFIXES:
            ad_type = line[77:79].strip().upper() if len(line) >= 78 else ''
            element = AD_TO_ELEMENT.get(ad_type, ad_type[:2].title() if ad_type else '')
            line = line[:66].ljust(76) + element.rjust(2)
        out.append(line)
    return out


def shift_negative_residues_lines(lines):
    """Shift residues per chain so all are >= 1.

    Two logical passes (need every chain's minimum before applying the shift),
    but each line is parsed only once: pass 1 records (index, chain, resnum)
    triples for atom lines while accumulating min_per_chain, pass 2 patches
    those indices in the output list.
    """
    parsed = []  # (line_index, chain, resnum) for atom lines
    min_per_chain = {}
    for idx, L in enumerate(lines):
        if L[:6] in ATOM_HETATM_PREFIXES and len(L) >= 26:
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

    # Copy and patch only the affected indices.
    out = list(lines)
    for idx, ch, rn in parsed:
        if ch in offsets:
            L = out[idx]
            out[idx] = L[:22] + f'{rn + offsets[ch]:>4}' + L[26:]
    return out


def extract_binding_site_lines(protein_lines, ligand_lines, radius=8.0):
    """Filter protein lines to only the residues with at least one atom
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
    for line in ligand_lines:
        if line[:6] in ATOM_HETATM_PREFIXES and len(line) >= 54:
            try:
                ligand_xyz.append((
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ))
            except ValueError:
                pass
    if not ligand_xyz:
        return protein_lines  # nothing to filter against

    r2 = radius * radius

    # First pass — identify (chain, resid) pairs to keep.
    keep = set()
    for line in protein_lines:
        if line[:6] not in ATOM_HETATM_PREFIXES or len(line) < 54:
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
        return protein_lines  # safety: keep everything if nothing matched

    # Second pass — output only matching residue lines plus essential headers.
    out = []
    for line in protein_lines:
        if line[:6] in ATOM_HETATM_PREFIXES and len(line) >= 26:
            chain = (line[21] if len(line) > 21 else '').strip() or 'A'
            resid = line[22:26].strip()
            if (chain, resid) in keep:
                out.append(line)
        elif line[:6] in PDB_HEADER_PREFIXES or line.startswith('TER') or line.startswith('END'):
            out.append(line)
    return out


def merge_protein_and_ligand_lines(protein_lines, ligand_lines):
    """Concat protein + ligand lines, with TER between, ending with END."""
    p = [L for L in protein_lines if not L.startswith('END')]
    l = [L for L in ligand_lines if L[:6] in MERGE_LIGAND_RECORDS]
    return p + ['TER'] + l + ['END']


def detect_ligand_resname_lines(lines):
    """Return the most-common non-water HETATM resname.

    Single pass; resname slice is computed once per line; lines too short to
    contain a resname are skipped explicitly rather than producing an empty
    one (PDB lines are normally 80 cols but some tools truncate trailing
    whitespace).
    """
    counts = Counter()
    for L in lines:
        if not L.startswith('HETATM') or len(L) < 20:
            continue
        rn = L[17:20].strip()
        if rn and rn not in SKIP_RESNAMES:
            counts[rn] += 1
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


# CSS + JS injected into ProLIF's vis.js HTML before returning to the caller:
#  - aggressive margin/padding reset on body and wrapper divs
#  - network container has a fixed height so the body's natural height is
#    exactly `network + legend` with no surplus
#  - after vis.js finishes its physics simulation, network.fit() zooms its
#    viewport to the actual graph extent (removes internal canvas whitespace),
#    then `parent.postMessage({type: 'prolif-ready', height: N})` reports the
#    exact content height — the TS side reads N and resizes the iframe to fit
_COMPACT_CSS_JS = """<style>
  * { box-sizing: border-box; }
  html { margin: 0 !important; padding: 0 !important; }
  body {
    margin: 0 !important; padding: 0 !important;
    background: white;
  }
  body > * { margin: 0 !important; }
  div[id^="mynetwork"], .vis-network, #mynetwork {
    height: 520px !important;
    width: 100% !important;
    padding: 0 !important;
    margin: 0 !important;
  }
  [class*="legend"], div[class*="-legend"] {
    padding: 4px 8px !important;
    margin: 0 !important;
  }
</style>
<script>
  document.addEventListener('DOMContentLoaded', function() {
    var sent = false;
    var send = function() {
      if (sent) return;
      sent = true;
      var h = Math.max(
        document.documentElement.scrollHeight || 0,
        document.body ? (document.body.scrollHeight || 0) : 0
      );
      parent.postMessage({ type: 'prolif-ready', height: h }, '*');
    };
    var tryHook = function() {
      if (typeof network !== 'undefined' && network && network.on) {
        network.on('stabilizationIterationsDone', function() {
          try { network.fit({animation: false}); } catch (e) {}
          setTimeout(send, 60);
        });
        setTimeout(send, 250);
      } else {
        setTimeout(tryHook, 50);
      }
    };
    tryHook();
  });
</script>"""


def make_html_compact(html: str) -> str:
    """Inject the compact CSS + vis.js readiness/resize hook into ProLIF's
    LigNetwork HTML so the result drops straight into a Datagrok iframe with
    no surrounding whitespace and self-reports its content height."""
    if '<head>' in html:
        return html.replace('<head>', '<head>' + _COMPACT_CSS_JS)
    if '<body>' in html:
        return html.replace('<body>', '<body>' + _COMPACT_CSS_JS)
    return _COMPACT_CSS_JS + html


# ---------- main ----------
# We splitlines() once per text input (here for `protein`, and below for
# `ligand` if/when it's available) and pass list[str] through every helper
# so the per-call cost of splitlines is paid exactly once per input.

# 1. Normalise protein input (PDBQT → PDB if needed; shift negative residues)
protein_lines = protein.splitlines()
if is_pdbqt_lines(protein_lines):
    protein_lines = pdbqt_to_pdb_lines(protein_lines)
protein_lines = shift_negative_residues_lines(protein_lines)

# 2. Detect ligand resname BEFORE merging — when a separate ligand PDB is
# provided (Docking case), trust IT for the resname rather than the receptor,
# whose HETATMs (buffers, terminal NH3+, ions) would otherwise outvote the
# real ligand via most-common heuristic.
ligand_provided = bool(ligand and ligand.strip())
ligand_lines = ligand.splitlines() if ligand_provided else []
if not ligand_resname:
    if ligand_provided:
        ligand_resname = detect_ligand_resname_lines(ligand_lines)
    else:
        ligand_resname = detect_ligand_resname_lines(protein_lines)

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
if ligand_provided and is_pdbqt_lines(ligand_lines):
    ligand_lines = pdbqt_to_pdb_lines(ligand_lines)

# 3b. BSV-case fix: if no separate ligand was provided (protein cell contains
# both protein and ligand), extract the matching HETATMs from the protein text
# now and remove them from the protein. This way pdbfixer only sees the
# protein, and we re-merge the original ligand HETATMs (with their original
# chain/resid intact) after pdbfixer.
if not ligand_provided:
    extracted_lines = []
    remaining_lines = []
    for line in protein_lines:
        # Need at least cols 1-26 for resname/chain/resid; reject shorter lines
        # so the slices below are always safe.
        if line.startswith('HETATM') and len(line) >= 26 \
                and line[17:20].strip() == sel_resname:
            ln_chain = line[21].strip() or 'A'
            ln_resid = line[22:26].strip()
            chain_ok = (sel_chain is None) or (ln_chain == sel_chain)
            resid_ok = (sel_resid is None) or (ln_resid == sel_resid)
            if chain_ok and resid_ok:
                extracted_lines.append(line)
                continue
        remaining_lines.append(line)
    if extracted_lines:
        ligand_lines = extracted_lines
        protein_lines = remaining_lines
        ligand_provided = True

# 4. pdbfixer for protein hydrogens — entirely in-memory (no temp files).
# Speed optimization: shrink the protein to just the binding-site neighborhood
# (8 Å around the ligand) before running pdbfixer. Hydrogenation scales
# linearly with atom count, so for a 300-residue protein this is typically a
# 5-10× speedup.
if ligand_provided:
    protein_for_fixer_lines = extract_binding_site_lines(
        protein_lines, ligand_lines, radius=8.0,
    )
else:
    protein_for_fixer_lines = protein_lines

# pdbfixer needs a text stream, so we join only at this boundary
protein_for_fixer = '\n'.join(protein_for_fixer_lines)
fixer = PDBFixer(pdbfile=StringIO(protein_for_fixer))
# We only need addMissingHydrogens() — the other pdbfixer steps don't apply:
#   findMissingResidues() — detects chain gaps; after binding-site extraction
#     we have intentional gaps and we DON'T want them filled. `missingResidues`
#     defaults to {} when skipped.
#   findMissingAtoms() / addMissingAtoms() — fixes missing side-chain atoms;
#     skipping saves ~100-500ms. Binding-site residues are usually
#     well-resolved in deposited structures.
fixer.addMissingHydrogens(pH=7.0)

# 5. Build the MDAnalysis Universe. The OpenMM-direct path
# `mda.Universe(fixer.topology, fixer.positions)` fails on some MDA versions
# because positions are treated as a filename — instead we serialize the
# H-fixed structure to PDB text in memory (no file I/O) and load that.
_buf = StringIO()
PDBFile.writeFile(fixer.topology, fixer.positions, _buf)
fixed_protein = _buf.getvalue()

# Merge ligand lines into the fixed protein (operating on lines avoids one
# more splitlines call — `merge_protein_and_ligand_lines` returns a list).
if ligand_provided:
    merged_lines = merge_protein_and_ligand_lines(
        fixed_protein.splitlines(), ligand_lines,
    )
    merged_pdb = '\n'.join(merged_lines)
else:
    merged_pdb = fixed_protein

u = mda.Universe(NamedStream(StringIO(merged_pdb), 'merged.pdb'))

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
    # `byres` selects the ENTIRE residue when any of its atoms is within the
    # cutoff, so atoms further out are still included.
    'protein and byres around 5.0 group ligand', ligand=ligand_ag
)
protein_ag.guess_bonds()
ligand_ag.guess_bonds()

n_h_lig = len(ligand_ag.select_atoms('element H'))

# 6. Prepare ligand mol — write to in-memory PDB via NamedStream, parse with
# RDKit. No temp file roundtrip.
if n_h_lig > 0:
    ligand_mol_obj = plf.Molecule.from_mda(ligand_ag)
else:
    _ligbuf = StringIO()
    ligand_ag.write(NamedStream(_ligbuf, 'lig.pdb'))
    lig_pdb = _ligbuf.getvalue()
    lig_rdkit = Chem.MolFromPDBBlock(
        lig_pdb, removeHs=False, sanitize=False, proximityBonding=False,
    )
    try:
        rdDetermineBonds.DetermineBonds(lig_rdkit, charge=0)
    except Exception:
        lig_rdkit = Chem.MolFromPDBBlock(
            lig_pdb, removeHs=False, sanitize=False, proximityBonding=True,
        )
    Chem.SanitizeMol(lig_rdkit)
    # CRITICAL: addResidueInfo=True propagates PDBResidueInfo to new H atoms —
    # without it ProLIF's residue splitting crashes silently on Windows.
    lig_rdkit = Chem.AddHs(lig_rdkit, addCoords=True, addResidueInfo=True)
    ligand_mol_obj = plf.Molecule.from_rdkit(lig_rdkit)

# Protein conversion. The MDA→RDKit converter occasionally hits
# AtomValenceException on residues whose bond perception produces an
# over-valent atom (modified residues, awkward pdbfixer geometry, glycans,
# etc.). Try progressively more permissive modes:
#   1. force=True   — skip strict valence check after bond inference
#   2. NoImplicit=False — skip bond inference entirely; use implicit Hs.
try:
    protein_mol_obj = plf.Molecule.from_mda(protein_ag, force=True)
except Exception as e_force:
    logger.warning('from_mda(force=True) failed: %s: %s',
                   type(e_force).__name__, e_force)
    logger.info('Retrying with NoImplicit=False (skip bond-order inference)')
    try:
        protein_mol_obj = plf.Molecule.from_mda(protein_ag, NoImplicit=False)
    except TypeError:
        protein_mol_obj = plf.Molecule.from_mda(protein_ag)

# 7. Fingerprint with Datagrok-aligned SMARTS where it makes sense.
# If our combined SMARTS fail to parse for any reason, fall back to ProLIF's
# defaults so the script still produces a useful result.
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
    logger.warning('Datagrok-aligned SMARTS rejected (%s); using ProLIF defaults', e)
    fp = plf.Fingerprint(interactions=interactions)
fp.run_from_iterable([ligand_mol_obj], protein_mol_obj, n_jobs=1)

# 8. LigNetwork → HTML. The TS layer injects the compact CSS + ready-hook on
# the way to the iframe, so we just return the raw ProLIF HTML here.
view = fp.plot_lignetwork(ligand_mol_obj, kind='frame', frame=0)
html = extract_html_from_view(view)
