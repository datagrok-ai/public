#top-menu: Chem | Calculate | IUPAC Name
#name: IUPAC Name
#description: Generates IUPAC names for molecules deterministically using openclatura (Blue Book 2013 rules).
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, pip, rdkit, {pip: [openclatura]}]
#meta.domain: chem
#input: dataframe table {nullable: false} [Input data table]
#input: column molecules {type: categorical; semType: Molecule; nullable: false}
#meta.timeout: 900000
#output: dataframe result {action:join(table)} [IUPAC names; empty where the molecule could not be named]

import pandas as pd
from rdkit import Chem
from openclatura import name_many

values = table[molecules].tolist()
indices = []
inputs = []
for i, value in enumerate(values):
    if not isinstance(value, str) or value == '':
        continue
    if 'M  END' in value:
        mol = Chem.MolFromMolBlock(value, sanitize=True)
        if mol is None:
            continue
        indices.append(i)
        inputs.append(mol)
    else:
        indices.append(i)
        inputs.append(value)

names = [''] * len(values)
if inputs:
    for i, r in zip(indices, name_many(inputs, processes='auto', verify_opsin=False)):
        if r.ok and r.name:
            names[i] = r.name

result = pd.DataFrame({'iupac_name': names})
