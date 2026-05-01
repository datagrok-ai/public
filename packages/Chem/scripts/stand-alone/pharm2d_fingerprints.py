#name: Pharm2D Fingerprints
#description: Computes Pharm2D (Gobbi-Poppinger 1998) topological pharmacophore fingerprints over the Chem package's standard 7-family SMARTS definitions, using RDKit's SigFactory with pair AND triplet feature combinations (maxPointCount=3) and triangle-pruned distance bins (0,2)(2,4)(4,6)(6,8)(8,12). Returns the sparse-bit-vector "on bits" as a comma-separated string per molecule for transport back to TypeScript.
#help-url: https://www.rdkit.org/docs/source/rdkit.Chem.Pharm2D.html
#language: python
#meta.domain: chem
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#input: dataframe features [Pharmacophore SMARTS definitions table — must have columns 'family' and 'smarts']
#output: dataframe fingerprints {action:join(data)} [Per-molecule Pharm2D fingerprint as comma-separated on-bit indices]

import inspect
import pandas as pd
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
# Alias to a unique name so this binding can't collide with anything cached
# from a previous script-runner pass (Datagrok reuses the Python process,
# and `SigFactory` is a Pharm2D submodule name as well as the class name —
# importing the bare name has been observed to resolve to the module on
# subsequent runs even when the source says `from .SigFactory import SigFactory`).
from rdkit.Chem.Pharm2D.SigFactory import SigFactory as _SigFactoryClass

if not inspect.isclass(_SigFactoryClass):
    raise TypeError(
        f'Expected Pharm2D.SigFactory.SigFactory to be a class, got '
        f'{type(_SigFactoryClass).__name__}: {_SigFactoryClass!r}')

# Family-letter → FDEF-friendly name. Mirrors FAMILY_LETTER_TO_ID in
# packages/Chem/src/panels/pharmacophore-features.ts:11-19. The letters come
# from pharmacophore-features.csv. FDEF feature names disallow spaces, so
# 'Halogen Bond' is renamed 'HalogenBond' here — purely an internal label,
# not visible to the user.
FAMILY_LETTER_TO_NAME = {
    'D': 'Donor',
    'A': 'Acceptor',
    'H': 'Hydrophobic',
    'a': 'Aromatic',
    'P': 'Positive',
    'N': 'Negative',
    'X': 'HalogenBond',
}

# Build the FDEF text from the supplied SMARTS table — one DefineFeature
# block per row; multiple SMARTS share the same Family if they belong to the
# same pharmacophore class (e.g. several Donor SMARTS).
#
# RDKit's FDEF parser requires the `Weights` line to have one weight per atom
# in the SMARTS pattern. Single-atom SMARTS get `1.0`, multi-atom ones get
# space-separated `1.0 1.0 1.0 ...` so the feature is placed at the atom
# centroid (the standard convention; alternatives like first-atom-only via
# `1 0 0 ...` are also valid). Pre-parsing each SMARTS via `MolFromSmarts`
# tells us the atom count and lets us silently skip any SMARTS RDKit can't
# parse rather than crashing the whole factory build.
def _smarts_atom_count(smarts):
    try:
        qmol = Chem.MolFromSmarts(smarts)
        return qmol.GetNumAtoms() if qmol is not None else 0
    except Exception:
        return 0


# Build the FDEF one feature at a time and validate each in isolation. The
# FDEF parser sometimes counts atoms differently from `MolFromSmarts` (e.g.
# for recursive SMARTS like `[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]`), so a
# single-feature `BuildFeatureFactoryFromString` is the only reliable check.
# We also use comma-separated weights, which matches RDKit's BaseFeatures.fdef
# convention (e.g. `Weights 1.0,1.0,1.0`).
def _try_build_one(feat_id, family_name, smarts, n_atoms):
    weights = ','.join(['1.0'] * n_atoms)
    block = (
        f'DefineFeature {feat_id} {smarts}\n'
        f'  Family {family_name}\n'
        f'  Weights {weights}\n'
        'EndFeature'
    )
    try:
        ChemicalFeatures.BuildFeatureFactoryFromString(block)
        return block, None
    except Exception as e:
        return None, str(e)


fdef_blocks = []
skipped = []
for i in range(len(features)):
    family_letter = features['family'].iloc[i]
    smarts = features['smarts'].iloc[i]
    if not family_letter or not smarts:
        continue
    family_name = FAMILY_LETTER_TO_NAME.get(family_letter)
    if not family_name:
        continue
    n_atoms = _smarts_atom_count(smarts)
    if n_atoms == 0:
        skipped.append(f'{family_letter}/{smarts} :: MolFromSmarts returned None')
        continue
    feat_id = f'{family_letter}_{i}'
    block, err = _try_build_one(feat_id, family_name, smarts, n_atoms)
    if block is None:
        # Try with single-weight fallback in case FDEF parser counts atoms
        # differently from MolFromSmarts (recursive-SMARTS edge cases).
        block2, err2 = _try_build_one(feat_id, family_name, smarts, 1)
        if block2 is not None:
            fdef_blocks.append(block2)
            continue
        skipped.append(f'{family_letter}/{smarts} :: {err}')
        continue
    fdef_blocks.append(block)

if skipped:
    print(f'Pharm2DFingerprints: {len(skipped)} broken features (skipped):')
    for s in skipped[:20]:
        print(f'  - {s}')
print(f'Pharm2DFingerprints: using {len(fdef_blocks)} working features out of {len(features)}')

fdef_lines = ['\n'.join(fdef_blocks)]

fdef_text = '\n'.join(fdef_lines)
feat_factory = ChemicalFeatures.BuildFeatureFactoryFromString(fdef_text)

# Pharm2D signature factory — pairs AND triplets (maxPointCount=3),
# triangle-pruned bins (drops geometrically impossible distance triples,
# keeping the bit count manageable). Bin edges follow the Gobbi-Poppinger
# 1998 reference and the RDKit Gobbi_Pharm2D default.
sig_factory = _SigFactoryClass(
    feat_factory,
    minPointCount=2,
    maxPointCount=3,
    trianglePruneBins=True,
)
sig_factory.SetBins([(0, 2), (2, 4), (4, 6), (6, 8), (8, 12)])
sig_factory.Init()


_first_error = [None]


def _smiles_to_onbits(smi):
    """Returns comma-separated on-bit indices of the Pharm2D fingerprint, or
    a diagnostic string starting with 'ERR:' on the first failure (so the
    caller can see what's actually breaking) and empty string on subsequent
    failures (to avoid spamming)."""
    if not smi:
        return ''
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return ''
        fp = Generate.Gen2DFingerprint(mol, sig_factory)
        on_bits = list(fp.GetOnBits())
        return ','.join(str(b) for b in on_bits)
    except Exception as e:
        if _first_error[0] is None:
            _first_error[0] = f'{type(e).__name__}: {e}'
            print(f'Pharm2DFingerprints: first per-molecule error on "{smi}": '
                  f'{_first_error[0]}')
            return f'ERR:{_first_error[0][:200]}'
        return ''


smiles_list = data[smiles].tolist()
fps = [_smiles_to_onbits(s) for s in smiles_list]

# If the first molecule produced an empty fingerprint, prepend a diagnostic
# string on the FIRST row so the caller can read it. This happens when the
# feature factory has no working features (all SMARTS rejected).
if fps and fps[0] == '':
    diag = (f'WORKING_FEATURES={len(fdef_blocks)} | '
            f'TOTAL_FEATURES={len(features)} | '
            f'SKIPPED={len(skipped)} | '
            f'FIRST_SKIPS={skipped[:3] if skipped else []} | '
            f'NUM_BITS={sig_factory.GetSigSize() if hasattr(sig_factory, "GetSigSize") else "?"}')
    fps[0] = f'DIAG:{diag}'

fingerprints = pd.DataFrame({'pharm2d_fp': fps})
