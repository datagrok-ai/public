#name: RGroupGetter
#language: python
#input: string smiles
#input: dataframe df1
#input: string core
#input: string prefix
#output: dataframe result

from rdkit import Chem
import numpy as np

result = {}
smiles = df1[smiles].tolist()
length = len(smiles)
core = Chem.MolFromSmiles(core)
if core is not None:
  fragments = dict()
  r_group = 1
  r_group_map = dict()
  for n in range(0, length):
    mol = Chem.MolFromSmiles(smiles[n])
    if mol is None:
      continue
    mol_no_core = Chem.ReplaceCore(mol, core, labelByIndex=True)
    if mol_no_core is not None:
      try:
        fragments_ = list(Chem.GetMolFrags(mol_no_core, asMols=True))
      except:
        fragments_ = []
    else:
      fragments_ = []
    for fragment in fragments_:
      fragment_smiles = Chem.MolToSmiles(fragment)
      id = fragment_smiles[:fragment_smiles.find(']')+1]
      if id not in r_group_map:
        r_group_map[id] = f'{prefix}{r_group}'
        fragments[r_group_map[id]] = np.full(length, None, dtype=object)
        r_group += 1
      fragments[r_group_map[id]][n] = fragment_smiles.replace(id, f'[{r_group_map[id]}]', 1)
else:
  raise Exception('Core is empty')
result = pd.DataFrame.from_dict(fragments)
