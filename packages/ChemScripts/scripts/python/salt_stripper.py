#name: Salt Stripper
#description: Removes salts from molecules and display the salt stripped molecules
#language: python
#sample: chem/smiles.csv
#tags: demo, chem, rdkit
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#output: dataframe stripped {action:join(data); semType: Molecule} [Stripped molecules, in SMILES format]

import numpy as np
from rdkit import Chem

smiles = data[smiles]

salts = [
    "[Cl,Br,I]",
    "[Li,Na,K,Ca,Mg]",
    "[O,N]",
    "[N](=O)(O)O",
    "[P](=O)(O)(O)O",
    "[P](F)(F)(F)(F)(F)F",
    "[S](=O)(=O)(O)O",
    "[CH3][S](=O)(=O)(O)",
    "c1cc([CH3])ccc1[S](=O)(=O)(O)",
    "FC(F)(F)C(=O)O",
    "OC(=O)C=CC(=O)O",
    "OC(=O)C(=O)O",
    "OC(=O)C(O)C(O)C(=O)O",
    "C1CCCCC1[NH]C1CCCCC1"
]
salts = [Chem.MolFromSmarts(salt) for salt in salts]

length = len(smiles)
stripped = np.full(length, None, dtype=object)
for n in range(0, length):
    mol = Chem.MolFromSmiles(smiles[n])
    if mol is None:
        continue
    mol_stripping = mol
    for salt in salts:
        if len(Chem.GetMolFrags(mol_stripping)) > 1:
            mol_stripped = Chem.DeleteSubstructs(mol_stripping, salt)
            if mol_stripped.GetNumAtoms() > 0:
                mol_stripping = mol_stripped
    stripped[n] = Chem.MolToSmiles(mol_stripping)

# Convert to Pandas DataFrame
stripped = pd.DataFrame(stripped, columns=['stripped'])
