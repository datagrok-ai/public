#name: Murcko Scaffolds
#description: Generation of Murcko scaffolds from a molecule
#help-url: https://datagrok.ai/help/domains/chem/functions/murcko-scaffolds
#language: python
#sample: chem/smiles.csv
#tags: demo, chem, rdkit
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#output: dataframe scaffolds {action:join(data); semType: Molecule} [Murcko scaffolds, in SMILES format]

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

smiles = data[smiles]

length = len(smiles)
scaffolds = np.full(length, None, dtype=object)
for n in range(0, length):
    mol = Chem.MolFromSmiles(smiles[n])
    if mol is None:
        continue
    scaffolds[n] = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)

# Convert to Pandas DataFrame
scaffolds = pd.DataFrame(scaffolds, columns=['scaffolds'])
