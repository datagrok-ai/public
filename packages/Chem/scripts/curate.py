#name: CurateChemStructures
#description: curating a molecules set for structural data homogenization
#top-menu: Chem | Curate...
#language: python
#sample: chem/chem_standards.csv
#tags: demo, chem, rdkit
#input: dataframe data [Input data table]
#input: column smiles  {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#input: bool kekulization = false
#input: bool normalization = false
#input: bool reionization = false
#input: bool neutralization = false
#input: bool tautomerization = false
#input: bool mainFragment = false
#output: dataframe curated {action:join(data); semType: Molecule} [Molecules, in SMILES format]

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

smiles = data[smiles]
length = len(smiles)

standardized = np.full(length, None, dtype=object)

def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

if tautomerization:
    enumerator = rdMolStandardize.TautomerEnumerator()

for n in range(0, length):
    mol = Chem.MolFromSmiles(smiles[n], sanitize = True)

    if mol is None or mol.GetNumAtoms() == 0:
        continue
    if tautomerization:
        mol = enumerator.Canonicalize(mol)
    if normalization:
        mol = rdMolStandardize.Normalize(mol)
    if reionization:
        mol = rdMolStandardize.Reionize(mol)
    if neutralization:
        neutralize_atoms(mol)
    if mainFragment:
        mol = rdMolStandardize.FragmentParent(mol)
    if kekulization:
        Chem.Kekulize(mol)

    standardized[n] = Chem.MolToSmiles(mol, kekuleSmiles = kekulization)

curated = pd.DataFrame(standardized, columns = ['curated_molecule'])


