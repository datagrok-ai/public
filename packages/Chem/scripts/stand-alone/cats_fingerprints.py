#name: CATS Fingerprints
#description: Computes CATS2D (Schneider et al. 1999) topological pharmacophore-pair descriptors over the Chem package's 7-family SMARTS definitions. CATS2D counts (familyA, familyB, topological_distance) triples and normalises by feature counts (Schneider scaling) to make the descriptor scale-invariant. Designed specifically for scaffold-hop retrieval at low Tanimoto similarity. Returns the 7×7×10 = 490-dim float vector per molecule, encoded as space-separated floats for transport.
#help-url: https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1521-3838(199910)18:5%3C424::AID-QSAR424%3E3.0.CO%3B2-Y
#language: python
#meta.domain: chem
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#input: dataframe features [Pharmacophore SMARTS definitions table — must have columns 'family' and 'smarts']
#output: dataframe fingerprints {action:join(data)} [Per-molecule CATS2D fingerprint as space-separated floats]

import numpy as np
import pandas as pd
from rdkit import Chem

# Family-letter → family name. Mirrors the mapping in
# packages/Chem/src/panels/pharmacophore-features.ts:11-19 — same as our
# Pharm2D script so the two descriptors share a feature definition.
FAMILY_LETTER_TO_NAME = {
    'D': 'Donor',
    'A': 'Acceptor',
    'H': 'Hydrophobic',
    'a': 'Aromatic',
    'P': 'Positive',
    'N': 'Negative',
    'X': 'HalogenBond',
}
FAMILIES = list(FAMILY_LETTER_TO_NAME.values())   # 7 families
N_FAMILIES = len(FAMILIES)
N_BINS = 10                                       # topological distances 0..9 bonds

# Pre-compile SMARTS query mols per family from the supplied features table.
# Multiple SMARTS per family are unioned (any match → atom belongs to that family).
family_qmols = {f: [] for f in FAMILIES}
for i in range(len(features)):
    family_letter = features['family'].iloc[i]
    smarts = features['smarts'].iloc[i]
    if not family_letter or not smarts:
        continue
    family_name = FAMILY_LETTER_TO_NAME.get(family_letter)
    if not family_name:
        continue
    qmol = Chem.MolFromSmarts(smarts)
    if qmol is not None:
        family_qmols[family_name].append(qmol)


def _cats2d_vector(smi):
    """Returns the CATS2D pharmacophore vector for a SMILES, normalised per
    Schneider 1999 (each (A, B, d) bin divided by count(A)+count(B) to make
    the descriptor scale-invariant). Empty list on parse / RDKit failure."""
    if not smi:
        return []
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return []
        n_atoms = mol.GetNumAtoms()
        if n_atoms == 0:
            return []

        # Identify pharm atoms per family by unioning all family-SMARTS matches.
        fam_atoms = {f: set() for f in FAMILIES}
        for fam, qmols in family_qmols.items():
            for q in qmols:
                for match in mol.GetSubstructMatches(q):
                    fam_atoms[fam].update(match)

        # Topological distance matrix (atom-atom shortest path in bonds).
        dist_mat = Chem.GetDistanceMatrix(mol)

        # Pairwise feature-pair × distance histogram, 7 × 7 × 10.
        # Vectorised over atom pairs: for each (familyA, familyB) pair we
        # slice the atom-atom distance matrix at the family's atom sets via
        # np.ix_, mask out self-pairs and out-of-range distances, and use
        # np.add.at (NOT `hist[..., d] += 1` — that would only count each
        # bin once when multiple atom-pairs collide on the same distance)
        # to scatter-accumulate into the 10-distance-bin slice. Replaces
        # the per-atom-pair Python loop, which was O(|A|·|B|) Python-level
        # iterations per family pair and dominated runtime on large
        # screening libraries (~120K iter/molecule worst case).
        hist = np.zeros((N_FAMILIES, N_FAMILIES, N_BINS), dtype=np.float32)
        for fa_idx, fa in enumerate(FAMILIES):
            atoms_a = fam_atoms[fa]
            if not atoms_a:
                continue
            a_arr = np.fromiter(atoms_a, dtype=np.int32)
            for fb_idx, fb in enumerate(FAMILIES):
                atoms_b = fam_atoms[fb]
                if not atoms_b:
                    continue
                b_arr = np.fromiter(atoms_b, dtype=np.int32)
                d_mat = dist_mat[np.ix_(a_arr, b_arr)].astype(np.int32)
                mask = (d_mat >= 0) & (d_mat < N_BINS) & (a_arr[:, None] != b_arr[None, :])
                np.add.at(hist[fa_idx, fb_idx], d_mat[mask], 1)

        # Schneider 1999 normalisation: divide each (A, B) histogram by
        # (count(A) + count(B)) so the descriptor doesn't grow with molecule
        # size. Empty families leave their slice at zero.
        fam_counts = [len(fam_atoms[f]) for f in FAMILIES]
        for fa_idx in range(N_FAMILIES):
            for fb_idx in range(N_FAMILIES):
                denom = fam_counts[fa_idx] + fam_counts[fb_idx]
                if denom > 0:
                    hist[fa_idx, fb_idx, :] /= denom

        return hist.flatten().tolist()
    except Exception:
        return []


smiles_list = data[smiles].tolist()
fps = []
for s in smiles_list:
    v = _cats2d_vector(s)
    fps.append(' '.join(f'{x:.4f}' for x in v) if v else '')

fingerprints = pd.DataFrame({'cats_fp': fps})
