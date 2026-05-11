#name: ProteinLigandInteractionDiagram
#description: Interactive 2D protein-ligand interaction diagram via ProLIF
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, rdkit, mdanalysis, prolif, pdbfixer, openmm, matplotlib, requests, pip, {pip: [openbabel-wheel]}]
#meta.cache: none
#input: string protein
#input: string ligand
#input: string ligand_resname
#output: dataframe result

import logging
import os
import tempfile
from collections import Counter
from io import StringIO
from itertools import islice

# Point urllib at certifi's CA bundle BEFORE importing pdbfixer (which uses
# urlopen at runtime to fetch CCD residue templates from RCSB). Datagrok's
# script-worker container has no system CA bundle, so without this the SSL
# handshake fails and pdbfixer silently swallows the error (bare `except:`)
# in `_downloadCCDDefinition` and skips template lookup entirely. The downstream
# effect is that addMissingHydrogens does NOT add Hs to non-standard ligands
# (E4Y, etc.) — the ligand keeps its 0-H state, the script falls into the
# RDKit fallback path, and DetermineBonds fails with a charge mismatch
# (the missing Hs make heavy atoms appear charged), producing all-SINGLE bonds.
# Setting SSL_CERT_FILE points urllib at certifi's bundle, the CCD download
# succeeds, pdbfixer adds Hs, and downstream bond perception works.
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

logger = logging.getLogger('prolif')


# Van der Waals radii dict used by `AtomGroup.guess_bonds()`. MDAnalysis's
# built-in table keys are UPPERCASE element symbols ("CL", "BR", "NA"),
# but PDB topologies (and pdbfixer output) routinely use mixed-case atom
# types ("Cl", "Br", "Na"). Without these aliases, guess_bonds raises:
#
#   ValueError: vdw radii for types: Cl. These can be defined manually
#   using the keyword 'vdwradii'
#
# Solution mirrors openfe_analysis's `prolif.py` fix
# (github.com/OpenFreeEnergy/openfe_analysis @ prolif_class_implementation):
# clone the MDA table and add the common mixed-case forms as aliases.
# Iodine ("I") + chlorine/bromine/sodium covers the cases we see in
# drug-like ligands and AutoDock poses.
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
# XBAcceptor also keeps ProLIF's default — `pharmacophore-features.csv`
# defines no Family-X acceptor pattern, so there's no Datagrok-side reference
# to align against (only the donor entry, feature 42, is in the CSV).
# ----------------------------------------------------------------------------

# Family D — HB donor (donor heavy atom + bonded H). ProLIF's HBDonor expects
# the SMARTS to match a 2-tuple (donor, H), so the trailing [H] is required.
# Mirrors `pharmacophore-features.csv` features 3-7 (Nitrogen, Oxygen, Sulfur,
# Charged N-H, Acid O-H donors).
DG_DONOR_SMARTS = (
    '[$([#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]),'
    '$([#8!H0&!$([OH][C,S,P]=O)]),'
    '$([#16!H0]),'
    '$([#7;+;!H0]),'
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
# any explicit + charge. Mirrors `pharmacophore-features.csv` features 14-20.
# `&!$(...)` exclusion guards prevent peptide-amide nitrogens and zwitterion
# atoms from being misclassified as cationic. The exclusion attaches to the
# OUTER atom expression — the recursive `$(...)` must close before `&`.
DG_POSITIVE_SMARTS = (
    '[$([NX3]([CX4])([CX4,#1])[CX4,#1])&!$([NX3]-*=[!#6]),'
    '$([CX3](=N)(-N)[!N]),'
    '$(N=[CX3](N)-N),'
    '$([+,+2,+3])&!$(*[-,-2,-3]),'
    '$([NX3;H2;+0][CX4])&!$([NX3;H2](C)C=[O,N,S]),'
    '$([NX3;H1;+0]([CX4])[CX4])&!$([NX3;H1](C)(C)C=[O,N,S]),'
    '$(c1c[nH]cn1)]'
)

# Family N — Anionic atom: tetrazole, sulfonate/phosphonate, carboxylate,
# sulfonamide, any explicit − charge. Mirrors `pharmacophore-features.csv`
# features 21-25. The `&!$(*[+,+2,+3])` guard on the bare anion entry
# prevents zwitterion atoms from being misclassified — recursive `$(...)`
# closes before `&` so the exclusion binds to the outer atom expression.
DG_NEGATIVE_SMARTS = (
    '[$(c1nn[nH1]n1),'
    '$([SX4,PX4](=O)(=O)[O-,OH]),'
    '$([CX3,SX3,PX3](=O)[O-,OH]),'
    '$([-,-2,-3])&!$(*[+,+2,+3]),'
    '$([NR;H1,-1](-C(=O))-[SX4](=O)(=O))]'
)

# Family X — Halogen-bond donor (carbon + bonded halogen). ProLIF's XBDonor
# expects the SMARTS to match a 2-tuple (heavy atom, halogen).
DG_HALOGEN_DONOR_SMARTS = '[#6][Cl,Br,I;X1]'


# Mapping ProLIF interaction names → short codes for the comma-separated
# per-row Interactions column (e.g. "GLY351_HD, LYS350_PI"). HBDonor and
# HBAcceptor get DIFFERENT codes (HD/HA) since direction matters chemically
# — a residue can be both donor and acceptor on different chemical groups,
# and the user may want to filter for one or the other.
INT_CODES = {
    'HBDonor': 'HD', 'HBAcceptor': 'HA',
    'Hydrophobic': 'HY',
    'PiStacking': 'PI',
    'Cationic': 'CAT', 'Anionic': 'AN',
    'CationPi': 'CPI', 'PiCation': 'CPI',
    'XBDonor': 'XB', 'XBAcceptor': 'XB',
    'MetalDonor': 'M', 'MetalAcceptor': 'M',
    'VdWContact': 'VDW',
}


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
    # waters
    'HOH', 'WAT', 'H2O', 'D2O', 'DOD',
    # ions / metals
    'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'MN', 'FE', 'CU', 'NI',
    # crystallographic buffers / additives
    'SO4', 'PO4', 'NO3', 'ACT', 'CO3',
    'GOL', 'EDO', 'PEG', 'PG4', 'DMS', 'TRS', 'IMD', 'BME',
    # common biological cofactors — typically not the inhibitor of interest
    # and they often outvote a small ligand on the most-common-HETATM
    # heuristic (cytochromes, dehydrogenases, kinases). Users can still
    # explicitly target a cofactor via `ligand_resname`.
    'HEM', 'HEC', 'HEB', 'HEA',                       # heme variants
    'NAD', 'NAI', 'NAP', 'NDP', 'NAH', 'NAJ',         # NAD(H)/NADP(H)
    'FAD', 'FMN', 'FDA',                              # flavins
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP',         # nucleotides
    'COA', 'ACO', 'COO',                              # CoA
    'SAM', 'SAH',                                     # SAM/SAH
    'PLP', 'PMP',                                     # PLP/PMP
    'BTN',                                            # biotin
    'CLA', 'CHL',                                     # chlorophyll
    'B12', 'COB', 'BCA',                              # cobalamin
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

# Records to preserve when concatenating ligand atoms into a fixed protein.
# We KEEP CONECT (the ligand's bonds are usually authoritative — from AutoDock
# or an upstream RDKit pass) but `merge_protein_and_ligand_lines` renumbers
# the ligand atom serials and updates the CONECT records so they don't
# collide with the protein's serials. Without renumbering, MDAnalysis would
# resolve a "CONECT 1 2" referring to the ligand to the protein's atoms 1,2
# instead, producing nonsense cross-bonds at the protein/ligand boundary.
MERGE_LIGAND_RECORDS = ('ATOM  ', 'HETATM', 'CONECT')


# Helper functions take/return list[str] (one element per PDB line) so the
# main flow only calls splitlines() once per text input. Calling splitlines on
# a few-thousand-line PDB is ~1-5ms each — small individually, but it adds up
# across 6-8 calls per row. Operating on lists everywhere makes the cost
# proportional to the work, not to the number of helpers.

def is_pdbqt_lines(lines):
    # Ligand PDBQT marker records — short-circuit on the FIRST ligand marker
    # seen before any ATOM/HETATM record. (Once we hit ATOM/HETATM without
    # having seen ROOT/BRANCH, this branch is decided.)
    for line in lines:
        if line.startswith(PDBQT_LIGAND_MARKERS):
            return True
        if line.startswith(ATOM_HETATM_PREFIXES):
            break
    # Receptor PDBQT: no ROOT/BRANCH, but ATOM/HETATM lines carry AutoDock atom
    # types in cols 78-79 (OA, NA, HD, NS, SA, HS, NX, OS) instead of standard
    # element symbols. Sample the first 50 ATOM/HETATM lines with `any()` over
    # an `islice`-capped generator — short-circuits on first match, no manual
    # counter, no None-padding.
    candidates = (line for line in lines
                  if line.startswith(ATOM_HETATM_PREFIXES) and len(line) >= 79)
    return any(line[77:79].strip().upper() in AUTODOCK_ATOM_TYPES
               for line in islice(candidates, 50))


def pdbqt_to_pdb_lines(lines):
    """Strip AutoDock-specific records and restore element symbols in cols 77-78."""
    out = []
    for line in lines:
        if line.startswith(PDBQT_RECORDS_TO_DROP):
            continue
        if line.startswith(ATOM_HETATM_PREFIXES):
            # Need len >= 79 to safely read the 2-char AutoDock atom-type
            # field at cols 78-79 (0-indexed [77:79]); a shorter line would
            # silently truncate 2-letter elements like 'CL' to 'C'.
            ad_type = line[77:79].strip().upper() if len(line) >= 79 else ''
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
    pdbfixer's findMissingResidues()/addMissingResidues() — otherwise
    pdbfixer would try to "fill" the gaps we deliberately removed with
    phantom amino acids. We DO call findMissingAtoms/addMissingAtoms
    on the kept residues themselves to fix missing side-chain atoms
    (without this, residues like a HIS missing one tautomer's atoms
    trigger "HIS residue has the wrong set of atoms" downstream in
    addMissingHydrogens).
    """
    # Collect ligand atom coordinates.
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
        return protein_lines  # nothing to filter against

    r2 = radius * radius

    # First pass — identify (chain, resid) pairs to keep.
    # Single sweep parses every atom line into parallel arrays; numpy then
    # computes the protein-vs-ligand min-distance matrix in one broadcast.
    # On a typical 3000-atom protein × 30-atom ligand this is ~3-30× faster
    # than the per-atom Python loop, with byte-identical output (verified
    # against scripts/bench_extract_binding_site.py).
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
        return protein_lines  # no atoms parsed

    prot_xyz = np.column_stack((xs, ys, zs))             # (N, 3)
    lig_xyz_np = np.asarray(ligand_xyz, dtype=float)     # (L, 3)

    # Squared distance from each protein atom to its NEAREST ligand atom.
    # Broadcast: (N, 1, 3) - (1, L, 3) -> (N, L, 3); sum-sq xyz; min over L.
    diff = prot_xyz[:, None, :] - lig_xyz_np[None, :, :]
    min_d2 = (diff * diff).sum(axis=-1).min(axis=1)      # (N,)

    in_range = np.flatnonzero(min_d2 < r2)
    keep = {(chains[i], resids[i]) for i in in_range}

    if not keep:
        return protein_lines  # safety: keep everything if nothing matched

    # Second pass — output only matching residue lines plus essential headers.
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
    """Concat protein + ligand lines, with TER between, ending with END.

    Renumbers the ligand atom serials (and any matching CONECT references)
    to start ABOVE the highest protein serial, so MDAnalysis can resolve
    every CONECT record unambiguously when it loads the merged PDB. Without
    this, the ligand's CONECT records (typically referencing atoms 1..M)
    would collide with the protein's atoms 1..M and produce wrong bonds.
    """
    # Pass 1 over protein: drop END lines, find the highest atom serial.
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

    # Pass 2 over ligand: build old->new serial remap from atom lines.
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

    # Pass 3 over ligand: rewrite atom serials and CONECT serial references.
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
            # PDB CONECT format: 6-char "CONECT", then up to 4 atom serials in
            # 5-char fields (cols 7-11, 12-16, 17-21, 22-26 — sometimes more).
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
                # Malformed CONECT — drop it rather than risk wrong bonds.
                continue
            l_kept.append(L)

    return p_kept + ['TER'] + l_kept + ['END']


def ensure_unique_atom_names(lines):
    """Rewrite atom names within each (chain, resid) so every atom has a
    unique name within its residue. AutoDock pose output uses generic names
    (`C`, `N`, `O`, ...) that COLLIDE within the same residue — pdbfixer's
    OpenMM PDB parser then silently drops all but the first occurrence under
    each (residue, atom_name) key, leaving a tiny fragment instead of the
    full ligand. The new names are formed as `<element><index>` (e.g.
    `C1, C2, ..., N1, N2, ...`), padded right to 4 chars to fit cols 13-16.
    Idempotent: re-running on already-unique names is harmless because each
    residue starts a fresh counter.
    """
    counters = {}
    out = []
    for L in lines:
        if not L.startswith(ATOM_HETATM_PREFIXES) or len(L) < 27:
            out.append(L)
            continue
        chain = L[21]
        resid = L[22:26]
        # Prefer the explicit element symbol (cols 77-78); fall back to the
        # current atom-name field if the element column is blank.
        elem = L[76:78].strip() if len(L) >= 78 else ''
        if not elem:
            elem = L[12:16].strip().rstrip('0123456789') or 'X'
        key = (chain, resid)
        n = counters.get(key, {}).get(elem, 0) + 1
        counters.setdefault(key, {})[elem] = n
        new_name = f'{elem}{n}'.ljust(4)[:4]
        out.append(L[:12] + new_name + L[16:])
    return out


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
#    then `parent.postMessage({type:'prolif-ready', pngDataUrl: N})` ships
#    the canvas back to the TS side, which stores it in `PL Diagram` for
#    the `rawPng` cell renderer to draw. Identical mechanism to ProLIF's
#    own `LigNetwork.save_png()` — the legend is a separate HTML block,
#    NOT part of the canvas, so it's naturally excluded from the PNG.
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
    // Post the canvas as a PNG data URL once vis.js has settled. The 250ms
    // delay after stabilizationIterationsDone lets the post-fit redraw
    // finish; without it the canvas can be captured mid-animation, producing
    // a half-drawn PNG (empirically 250ms is enough on 30-residue networks).
    var send = function() {
      if (sent) return;
      sent = true;
      var pngDataUrl = null;
      try {
        var canvas = document.getElementsByTagName('canvas')[0];
        if (canvas && canvas.width > 0 && canvas.height > 0)
          pngDataUrl = canvas.toDataURL('image/png');
      } catch (e) { /* canvas missing — leave pngDataUrl null */ }
      try { parent.postMessage({ type: 'prolif-ready', pngDataUrl: pngDataUrl }, '*'); }
      catch (e) { /* parent gone (iframe detached) */ }
    };
    var tryHook = function() {
      if (typeof network !== 'undefined' && network && network.on) {
        network.on('stabilizationIterationsDone', function() {
          try { network.fit({animation: false}); } catch (e) {}
          setTimeout(send, 250);
        });
        // safety fallback if stabilization event never fires (some
        // single-residue networks skip physics simulation)
        setTimeout(send, 2000);
      } else {
        setTimeout(tryHook, 50);
      }
    };
    tryHook();
  });
</script>"""


def make_html_compact(html: str) -> str:
    """Inject the compact CSS + vis.js readiness hook into ProLIF's LigNetwork
    HTML so the result drops straight into a Datagrok iframe with no
    surrounding whitespace, and postMessages the canvas PNG once vis.js
    settles (consumed by the TS-side `_captureOne`)."""
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

# Make ligand atom names unique within each residue. AutoDock-style poses
# use generic names (every C is just "C", every N is just "N") which collide
# under pdbfixer's (residue, atom_name) deduplication and silently lose
# all but the first occurrence per element — the script then sees a 6-atom
# stub (1 C, 1 N, 1 O + some Hs) instead of the full ligand. See
# `ensure_unique_atom_names` for details.
if ligand_provided:
    ligand_lines = ensure_unique_atom_names(ligand_lines)

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

# 4. Merge protein + ligand BEFORE pdbfixer so addMissingHydrogens sees the
# WHOLE merged structure and adds Hs to BOTH the protein (template-based)
# and the ligand (template-based for standard PDB residues like HEM,
# geometry-based via OpenBabel in the no-Hs fallback below).
#
# Speed optimization: shrink the PROTEIN to the binding-site neighborhood
# (8 Å around the ligand) before merging. Hydrogenation scales linearly
# with atom count, so on a 300-residue protein this is typically a 5-10×
# speedup on addMissingHydrogens. (We previously suspected this fragmented
# the topology and confused pdbfixer's CCD lookup for non-standard ligands,
# but the real culprit was duplicate atom names in AutoDock poses — fixed
# by `ensure_unique_atom_names` above. With unique names, binding-site
# extraction is safe.)
if ligand_provided:
    protein_for_fixer_lines = extract_binding_site_lines(
        protein_lines, ligand_lines, radius=8.0,
    )
    structure_for_fixer_lines = merge_protein_and_ligand_lines(
        protein_for_fixer_lines, ligand_lines,
    )
else:
    structure_for_fixer_lines = protein_lines

# pdbfixer needs a text stream, so we join only at this boundary
structure_for_fixer = '\n'.join(structure_for_fixer_lines)
fixer = PDBFixer(pdbfile=StringIO(structure_for_fixer))
# Standardise residues + fill in missing side-chain atoms BEFORE
# addMissingHydrogens. Without this, residues with unusual atom sets
# (e.g. a HIS with one tautomer's atoms missing, or a non-standard
# variant like HSD/HSE/HSP labelled as HIS) trigger:
#
#   ValueError: HIS residue (21) has the wrong set of atoms
#
# from addMissingHydrogens, because OpenMM's template matcher can't
# resolve the residue. The four-step sequence below is the canonical
# pdbfixer recipe for normalising a deposited PDB:
#   1. findNonstandardResidues  — flag e.g. HSD, MSE, SEC...
#   2. replaceNonstandardResidues — rewrite to the standard parent
#   3. findMissingAtoms        — locate side-chain atoms expected but
#      absent (the actual fix for the HIS-wrong-atom-set error)
#   4. addMissingAtoms         — build them from the residue template
#
# We deliberately STILL skip findMissingResidues/addMissingResidues —
# after binding-site extraction we have intentional chain gaps and
# don't want pdbfixer filling them with phantom amino acids. But
# `findMissingAtoms()` reads `self.missingResidues` internally to
# decide which residues to inspect, so we have to set it explicitly
# to an empty dict — otherwise it raises `AttributeError: 'PDBFixer'
# object has no attribute 'missingResidues'`.
#
# Each cleanup step (non-standard residue replacement, missing-atom
# repair) is wrapped in its own try/except: pdbfixer occasionally
# raises low-level errors like `KeyError: 4` on unusual atom indices
# or `IndexError` on truncated residue lists, especially for
# crystallographic outliers (alt-locs, terminal cleavage, partial
# sequences). When that happens we'd rather skip the cleanup step
# than fail the whole row — `addMissingHydrogens` alone is enough to
# get a usable structure for the majority of PDBs. The added
# ~100-500ms cost is fine for correctness on the rows where it works.
fixer.missingResidues = {}
try:
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        fixer.replaceNonstandardResidues()
except Exception as _e_nsr:
    logger.warning('pdbfixer replaceNonstandardResidues failed (%s: %s); '
                   'continuing without residue standardisation',
                   type(_e_nsr).__name__, _e_nsr)
try:
    fixer.findMissingAtoms()
    if fixer.missingAtoms or fixer.missingTerminals:
        fixer.addMissingAtoms()
except Exception as _e_ma:
    logger.warning('pdbfixer addMissingAtoms failed (%s: %s); '
                   'continuing without missing-atom repair',
                   type(_e_ma).__name__, _e_ma)
fixer.addMissingHydrogens(pH=7.0)

# 5. Serialize the H-fixed merged structure to PDB text in memory and load it.
# PDBFile.writeFile renumbers atoms in topology order starting at 1, and
# emits CONECT records for the ligand's non-standard residue bonds. The
# protein's amino-acid bonds aren't emitted (PDBFile.writeFile expects
# readers to reconstruct them from residue templates) — MDAnalysis doesn't
# carry templates, but `protein_ag.guess_bonds()` below regenerates them
# from element distances, so we don't need explicit CONECT for the protein.
_buf = StringIO()
PDBFile.writeFile(fixer.topology, fixer.positions, _buf)
merged_pdb = _buf.getvalue()

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
    # 6.5 Å covers ProLIF's longest interaction cutoff: PiStacking is an
    # umbrella over FaceToFace (5.5 Å) AND EdgeToFace (6.5 Å), the latter
    # being the true longest. `byres` selects the ENTIRE residue when any
    # of its atoms is within the cutoff, so atoms further out are still
    # included. (ProLIF's own default `vicinity_cutoff` is 6.0 Å; we use
    # 6.5 to leave no edge-to-face stacking interactions on the table.)
    'protein and byres around 6.5 group ligand', ligand=ligand_ag
)
# Bond loading. Always re-guess bonds from distances — matches the reference
# validation script. CONECT-derived bonds (when present) can be incomplete
# or stale; running `guess_bonds()` regenerates clean connectivity from
# interatomic distances + element radii, which is what MDA's RDKitConverter
# needs to perceive bond orders correctly downstream.
#
# Pass `vdwradii=_VDW_RADII` (mixed-case-aliased MDA table) so ligands with
# halogens or metal ions don't blow up with "vdw radii for types: Cl" —
# see the comment on `_VDW_RADII` at module scope for the full story.
protein_ag.guess_bonds(vdwradii=_VDW_RADII)
ligand_ag.guess_bonds(vdwradii=_VDW_RADII)

n_h_lig = len(ligand_ag.select_atoms('element H'))

# 6. Prepare ligand mol. Two paths matching the reference validation script:
#   * Hs present (the normal case after pdbfixer ran): use MDA's converter
#     directly. The Hs encode valences which MDA's RDKitConverter uses to
#     assign bond orders correctly.
#   * No Hs (rare — ligand bypassed pdbfixer): write to tempfile, parse
#     with RDKit, run xyz2mol-style bond perception via `DetermineBonds`,
#     then sanitize + AddHs. `addResidueInfo=True` is critical here:
#     ProLIF's residue splitting silently breaks on Windows without it.
ligand_mol_obj = None
if n_h_lig > 0:
    try:
        ligand_mol_obj = plf.Molecule.from_mda(ligand_ag)
    except Exception as _e_mda:
        # MDA's RDKitConverter sometimes produces invalid Lewis structures
        # for PDB ligands with explicit Hs — typically when the input PDB has
        # an H atom bonded to two heavy atoms (CCD-template / pdbfixer
        # mismatch on residues like NAD with mixed protonation states).
        # Fall through to the OpenBabel path which rebuilds bonds + Hs from
        # coordinates alone, sidestepping the problematic input Hs entirely.
        logger.warning('plf.Molecule.from_mda(ligand_ag) failed (%s: %s); '
                       'falling through to OpenBabel path',
                       type(_e_mda).__name__, _e_mda)
        # Strip Hs from ligand_ag so OpenBabel rebuilds them from scratch
        ligand_ag = ligand_ag.select_atoms('not element H')
if ligand_mol_obj is None:
    # No-Hs path: AutoDock pose, any non-standard ligand pdbfixer can't
    # template-match, or fallback after MDA's RDKitConverter rejected the
    # explicit-Hs input. We need bond orders + Hs from coordinates alone.
    # RDKit's `rdDetermineBonds.DetermineBonds` can do this in principle
    # (xyz2mol algorithm), but it FAILS when the apparent valences imply a
    # large net charge that doesn't match the requested `charge` parameter
    # — which happens for many docking poses (carbonyls/amides look like
    # charged O/N atoms because they're missing their Hs). Open Babel's
    # `PerceiveBondOrders` works on the same input because it uses bond
    # length heuristics rather than charge balancing, and produces correct
    # double bonds + aromatic rings on standard drug-like ligands.
    # NamedTemporaryFile + delete=False is the canonical replacement for
    # mkstemp + os.close + os.unlink — same semantics (we own the file,
    # close it ourselves so OpenBabel can read it), no manual fd dance.
    # We can't use delete=True with Windows: OpenBabel re-opens the path
    # by name, which fails while the original handle is still open on
    # Windows. delete=False + explicit unlink in the finally block works
    # cross-platform.
    _lig_pdb_text = ''
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as _tmp:
        _ligpath = _tmp.name
    try:
        ligand_ag.write(_ligpath)
        with open(_ligpath, 'r') as _f:
            _lig_pdb_text = _f.read()
        from openbabel import pybel as _pybel
        _ob_mol = next(_pybel.readfile('pdb', _ligpath))
        _ob_mol.OBMol.PerceiveBondOrders()
        _ob_mol.addh()
        # Hand off to RDKit via MOL block (V2000) — preserves bond orders
        # and 3D coordinates that ProLIF's `from_rdkit` then consumes.
        _mol_block = _ob_mol.write('mol')
    finally:
        try:
            os.unlink(_ligpath)
        except OSError:
            pass
    # Multi-tier RDKit parse with increasing leniency. Each tier wraps the
    # parse + sanitize so we don't lose the row on errors like
    #
    #   AtomValenceException: Explicit valence for atom # 2 C, 5, is
    #   greater than permitted
    #
    # which OpenBabel occasionally triggers when it over-perceives a bond
    # order (typically an aromatic / amide carbon ending up with 5 bonds).
    # Tiers run from "give me a correct Lewis structure" to "give me
    # ANYTHING I can pass to ProLIF" — the further down the cascade we
    # fall, the more interactions degrade to geometry-only detection
    # (positions are always preserved; bond-order-dependent SMARTS like
    # HBDonor/HBAcceptor may miss), but Hydrophobic/VdW/PiStacking still
    # register fine.
    lig_rdkit = None
    # Tier 1 — OpenBabel MOL block, full sanitisation. Most rows.
    try:
        lig_rdkit = Chem.MolFromMolBlock(_mol_block, removeHs=False, sanitize=True)
    except Exception as _e:
        logger.warning('MolFromMolBlock (strict) failed: %s: %s', type(_e).__name__, _e)
    # Tier 2 — OpenBabel MOL block, lenient sanitisation (skip the
    # PROPERTIES bit so over-bonded carbons are tolerated). Same bonds
    # OpenBabel proposed, just without RDKit re-checking maximum valence.
    if lig_rdkit is None:
        try:
            _m = Chem.MolFromMolBlock(_mol_block, removeHs=False, sanitize=False)
            if _m is not None:
                Chem.SanitizeMol(_m, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
                lig_rdkit = _m
        except Exception as _e:
            logger.warning('MolFromMolBlock (lenient) failed: %s: %s', type(_e).__name__, _e)
    # Tier 3 — PDB proximity bonding. All single bonds (RDKit's
    # `proximityBonding` doesn't infer bond orders), but the molecule is
    # structurally intact and geometry-based interactions still register.
    # Try strict sanitisation first, then lenient.
    if lig_rdkit is None:
        try:
            _m = Chem.MolFromPDBBlock(
                _lig_pdb_text, removeHs=False, sanitize=False, proximityBonding=True,
            )
            if _m is not None:
                try:
                    Chem.SanitizeMol(_m)
                except Exception:
                    Chem.SanitizeMol(_m, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
                lig_rdkit = Chem.AddHs(_m, addCoords=True, addResidueInfo=True)
        except Exception as _e:
            logger.warning('MolFromPDBBlock fallback failed: %s: %s', type(_e).__name__, _e)
    if lig_rdkit is None:
        raise RuntimeError(
            'Could not build an RDKit ligand mol from any parsing path '
            '(OpenBabel MOL strict/lenient, PDB proximity-bonding). The '
            'ligand topology may be too unusual — check the input PDB.'
        )
    ligand_mol_obj = plf.Molecule.from_rdkit(lig_rdkit)

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

# 8a. LigNetwork → HTML for the context panel (interactive vis.js diagram).
# Inject the compact CSS + ready-hook here so the panel widget gets
# consistent styling without having to re-apply it on the TS side. The
# injected script also postMessage's the canvas-rendered PNG back to the
# TS side as soon as vis.js stabilizes — that's what populates the
# per-cell PNG in the `PL Diagram` grid column.
view = fp.plot_lignetwork(ligand_mol_obj, kind='frame', frame=0)
html = make_html_compact(extract_html_from_view(view))

# 9. Extract structured per-residue interactions from the Fingerprint so the
# TS layer can populate sortable/filterable metric columns alongside the
# diagram. Format chosen for grid filters: comma-separated `RESNAME+RESID_CODE`
# pairs (e.g. "GLY351_H, LYS350_PI"). User can filter for "GLY351" to find
# rows engaging that residue, or for "_PI" to find rows with pi-stacking, etc.
# fp.to_dataframe() returns a DataFrame indexed by Frame with a MultiIndex
# (ligand, protein_residue, interaction) on the columns. Each cell is bool.
# `fp.to_dataframe()` raises a pandas KeyError for the 3-level MultiIndex
# names ('ligand', 'protein', 'interaction') when no interactions were
# detected (the resulting DataFrame is empty and has no MultiIndex). Treat
# that as zero interactions instead of letting it abort the whole script.
unique_pairs = set()  # (residue_str, code) — dedups HBDonor + HBAcceptor on
                      # the same residue into a single 'H' contribution
type_counts = Counter()
try:
    fp_df = fp.to_dataframe()
    cols_iter = fp_df.columns if fp_df is not None else []
except Exception as _e_fpdf:
    logger.warning('fp.to_dataframe() failed (%s: %s); reporting zero interactions',
                   type(_e_fpdf).__name__, _e_fpdf)
    cols_iter = []
for col in cols_iter:
    try:
        if not bool(fp_df[col].iloc[0]):
            continue
        # col is a 3-tuple (ligand_id, protein_residue, interaction_type)
        _, residue_str, int_type = col
    except Exception:
        continue
    # ProLIF residue format is "RESNAME+RESID.CHAIN" (e.g. "TYR23.A").
    # Drop the chain suffix — the user-facing column omits it for compactness.
    resname_resid = str(residue_str).split('.')[0]
    code = INT_CODES.get(int_type, str(int_type)[:3].upper())
    pair = (resname_resid, code)
    if pair in unique_pairs:
        continue
    unique_pairs.add(pair)
    type_counts[code] += 1

interactions_str = ', '.join(
    f'{r}_{c}' for r, c in sorted(unique_pairs)
)

# Per-code counts. Each `type_counts[code]` is "number of distinct residues
# with at least one interaction of this code". `unique_pairs` already de-
# duplicates (residue, code) pairs, so a single residue contributes once
# per code. HBDonor and HBAcceptor get separate codes (HD/HA), so a residue
# that's both a donor and an acceptor (e.g. SER hydroxyl on different
# chemical groups) contributes 1 to each. The TS layer uses these output
# column names verbatim; keep them stable.
RESULT_CODE_COLUMNS = (
    ('n_hbond_donor',    'HD'),
    ('n_hbond_acceptor', 'HA'),
    ('n_hydrophobic',    'HY'),
    ('n_pistacking',     'PI'),
    ('n_cationic',       'CAT'),
    ('n_anionic',        'AN'),
    ('n_cationpi',       'CPI'),
    ('n_xbond',          'XB'),
    ('n_metal',          'M'),
    ('n_vdw',            'VDW'),
)

# Pack into a 1-row DataFrame so Datagrok routes each column back to the TS
# caller. The panel reads `result.col('html')`; the batch handler reads the
# metric columns to populate filterable dataframe columns. One column per
# ProLIF interaction family. HBDonor and HBAcceptor are kept separate
# (n_hbond_donor / n_hbond_acceptor); the remaining symmetric pairs
# (CationPi/PiCation, XBDonor/XBAcceptor, MetalDonor/MetalAcceptor) still
# collapse via INT_CODES so a residue is counted once per family.
result_columns = {
    'html': [html],          # interactive vis.js LigNetwork; the panel
                             # mounts it directly with the legend visible.
                             # Browser-side capture (in `prolif-panel.ts`)
                             # extracts the canvas PNG via postMessage and
                             # writes it to the `PL Diagram` column for the
                             # `rawPng` cell renderer to draw.
    'interactions': [interactions_str],
}
for _col_name, _code in RESULT_CODE_COLUMNS:
    result_columns[_col_name] = [type_counts.get(_code, 0)]
result_columns['n_total'] = [sum(type_counts.values())]
result = pd.DataFrame(result_columns)
