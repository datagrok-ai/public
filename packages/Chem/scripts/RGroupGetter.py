#name: RGroupGetter
#language: python
#input: string smiles
#input: dataframe df1
#input: string core
#input: string prefix
#output: dataframe result

from rdkit import Chem
import numpy as np
import re

def np_none(shape):
  return np.full(shape, None, dtype=object)

result = {}
smiles = df1[smiles].tolist()
length = len(smiles)
core = Chem.MolFromSmiles(core)
if core is not None:
  fragments = np_none(length)
  max_fragments_length = 0
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
    fragments[n] = [Chem.MolToSmiles(fragment) for fragment in fragments_]
    for m in range(0, len(fragments[n])):
      fragments[n][m] = re.sub(r"\[\d+\*\]", '[R%d]' % (m + 1), fragments[n][m])
      fragments[n][m] = re.sub(r"\*", '[R%d]' % (m + 1), fragments[n][m])
    if len(fragments[n]) > max_fragments_length:
      max_fragments_length = len(fragments[n])
  r_groups = np_none(max_fragments_length)
  for r in range(0, max_fragments_length):
    r_groups[r] = np_none(length)
  for n in range(0, length):
    fragments_ = fragments[n]
    if fragments_ is not None:
      for g in range(0, len(fragments_)):
        r_groups[g][n] = fragments_[g]
  for g in range(0, len(r_groups)):
    result[prefix + str(g + 1)] = r_groups[g].tolist()
else:
  raise Exception('Core is empty')
result = pd.DataFrame.from_dict(result)