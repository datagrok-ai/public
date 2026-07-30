#name: Curate
#description: Standardizes the dataset
#top-menu: Chem | Transform | Curate...
#language: python
#sample: chem/chem_standards.csv
#meta.domain: chem
#input: dataframe data {nullable: false} [Input data table]
#input: column molecules {type:categorical; semType: Molecule; nullable: false}
#input: bool kekulization = false [Write the result as a Kekule SMILES instead of an aromatic one]
#input: bool normalization = true [Apply RDKit functional-group normalization]
#input: bool reionization = true [Move protons to the strongest acidic sites]
#input: bool neutralization = true [Neutralize formal charges where possible]
#input: bool tautomerization = false [Canonicalize the tautomer — slow, and opinionated]
#input: bool mainFragment = true [Keep only the largest fragment, dropping salts and solvents]
#output: dataframe curated {action:join(data); semType: Molecule} [Molecules, in SMILES format]

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

molecules = data[molecules]
length = len(molecules)

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

    if molecules[n] == "" or molecules[n] != molecules[n]:
        standardized[n] = ""
        continue

    mol = Chem.MolFromMolBlock(molecules[n], sanitize = True) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n], sanitize = True)

    if mol is None or mol.GetNumAtoms() == 0:
        standardized[n] = ""
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
