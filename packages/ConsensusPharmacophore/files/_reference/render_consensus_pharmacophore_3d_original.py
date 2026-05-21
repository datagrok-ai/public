#name: Render Consensus Pharmacophore (3D)
#description: Renders a consensus pharmacophore model as a 3D structure for the Biostructure viewer. Each feature family becomes a colored pseudo-atom; B-factor encodes frequency.
#language: python
#tags: chem, demo
#input: dataframe consensus_model {caption: Consensus model} [Output of Ligand-Based 3D Pharmacophore: family, x, y, z, frequency, n_ligands]
#input: string ref_pdb_id = "none" { caption: Reference PDB ID; optional: true } [Leave as 'none' to skip; enter a PDB ID to overlay it (only aligns if it shares the consensus coordinate frame)]
#output: string pharmacophore3d {semType: Molecule3D} [Consensus pharmacophore as a PDB block]
#environment: channels: [conda-forge], dependencies: [python=3.11, pandas, requests]
#meta.queueName: python_docker

# ---------------------------------------------------------------------------
# REFERENCE ONLY — do not register. Kept for design provenance.
# Stage 5b is implemented in TypeScript at src/renderer.ts; this Python file
# is the original implementation that the TS port was derived from.
# Blueprint reference: Q17, Appendix A.
# ---------------------------------------------------------------------------

import pandas as pd
import requests

# Map pharmacophore family -> (PDB element, 3-letter resName).
# Elements are chosen so Mol*'s default CPK coloring already separates the
# families before any manual styling -- good for a cold-open demo:
#   Donor -> N (blue)   Acceptor -> O (red)   Aromatic -> S (yellow)
#   Hydrophobe -> C (grey)   PosIonizable -> NA (violet)   NegIonizable -> CL (green)
FAMILY_MAP = {
    'Donor':            ('N',  'HBD'),
    'Acceptor':         ('O',  'HBA'),
    'Aromatic':         ('S',  'ARO'),
    'Hydrophobe':       ('C',  'HYD'),
    'LumpedHydrophobe': ('C',  'HYD'),
    'PosIonizable':     ('NA', 'POS'),
    'NegIonizable':     ('CL', 'NEG'),
}
DEFAULT_FEATURE = ('FE', 'PH4')  # fallback for any unmapped family


def resolve_family(fam):
    key = str(fam).strip()
    if key in FAMILY_MAP:
        return FAMILY_MAP[key]
    low = key.lower()
    for k, v in FAMILY_MAP.items():
        if k.lower() == low:
            return v
    if 'donor' in low:    return FAMILY_MAP['Donor']
    if 'accept' in low:   return FAMILY_MAP['Acceptor']
    if 'arom' in low:     return FAMILY_MAP['Aromatic']
    if 'hydroph' in low:  return FAMILY_MAP['Hydrophobe']
    if 'pos' in low:      return FAMILY_MAP['PosIonizable']
    if 'neg' in low:      return FAMILY_MAP['NegIonizable']
    return DEFAULT_FEATURE


def hetatm(serial, elem, res, chain, res_seq, x, y, z, occ, bfac):
    """One fixed-width PDB HETATM record (columns per the PDB v3.3 spec)."""
    atom_name = elem if len(elem) == 2 else ' ' + elem  # cols 13-16
    return (
        f"HETATM{serial:>5} {atom_name:<4} {res:<3} {chain}{res_seq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}          {elem:>2}"
    )


freq_col = 'frequency' if 'frequency' in consensus_model.columns else None

ph4_lines = []
serial = 1
for _, row in consensus_model.reset_index(drop=True).iterrows():
    elem, res = resolve_family(row['family'])
    bfac = float(row[freq_col]) if freq_col else 1.0
    ph4_lines.append(hetatm(
        serial, elem, res, 'P', serial,
        float(row['x']), float(row['y']), float(row['z']),
        1.00, bfac,
    ))
    serial += 1
ph4_block = "\n".join(ph4_lines)

# Optional reference structure. Caveat: this only lines up with the
# pharmacophore if the upstream pipeline superposed structures onto this
# reference's native frame. Any fetch error falls back to pharmacophore-only.
ref_block = ""
ref_pdb_id = (ref_pdb_id or "").strip().upper()
if ref_pdb_id and ref_pdb_id not in {"NONE", "NA", "NULL", "-", ""}:
    try:
        resp = requests.get(
            f"https://files.rcsb.org/download/{ref_pdb_id}.pdb", timeout=30)
        resp.raise_for_status()
        keep = ("ATOM", "HETATM", "TER", "SSBOND", "HELIX", "SHEET", "CONECT")
        ref_block = "\n".join(
            ln for ln in resp.text.splitlines() if ln.startswith(keep)
        ) + "\n"
    except Exception:
        ref_block = ""

pharmacophore3d = (
    "REMARK   Consensus pharmacophore -- chain P, B-factor = frequency\n"
    + ref_block
    + ph4_block
    + "\nEND\n"
)
