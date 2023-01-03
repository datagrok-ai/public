#name: FindMCS
#language: python
#input: string molecules
#input: dataframe df
#input: bool returnSmarts
#output: string result

from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np

def np_none(shape):
  return np.full(shape, None, dtype=object)

molecules = df[molecules].tolist()
length = len(molecules)
mols = np_none(length)
idx_err = []
for n in range(0, length):
  try:
  	mol = Chem.MolFromMolBlock(molecules[n], sanitize = True) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n], sanitize = True)
  except:
    mol = None
  if mol is None:
    idx_err.append(n)
    continue
  mols[n] = mol
mols = mols[mols != np.array(None)]
result = rdFMCS.FindMCS(mols).smartsString if returnSmarts else Chem.MolToSmiles(Chem.MolFromSmarts(rdFMCS.FindMCS(mols).smartsString))