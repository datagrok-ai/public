#name: Desc
#language: python
#input: string smiles
#input: dataframe df1
#input: string selected
#input: dataframe df2
#output: dataframe out

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import pandas as pd
import multiprocessing

smiles = df1[smiles]
descriptors = df2[selected]

def np_none(shape):
  return np.full(shape, None, dtype=object)

def _get_descriptors_list():
  descs = Descriptors._descList[:]
  descs_extra = [
    ("PMI1", Chem.Descriptors3D.PMI1),
    ("PMI2", Chem.Descriptors3D.PMI2),
    ("PMI3", Chem.Descriptors3D.PMI3),
    ("NPR1", Chem.Descriptors3D.NPR1),
    ("NPR2", Chem.Descriptors3D.NPR2),
    ("RadiusOfGyration", Chem.Descriptors3D.RadiusOfGyration),
    ("InertialShapeFactor", Chem.Descriptors3D.InertialShapeFactor),
    ("Eccentricity", Chem.Descriptors3D.Eccentricity),
    ("Asphericity", Chem.Descriptors3D.Asphericity),
    ("SpherocityIndex", Chem.Descriptors3D.SpherocityIndex),
   ]
  for de in descs_extra:
    descs.append(de)
  return descs

def _get_descriptors_funcs(descriptors):
  descriptors_funcs = {}
  descs = _get_descriptors_list()
  for desc in descs:
    if desc[0] in descriptors:
      descriptors_funcs[desc[0]] = desc[1]
  return descriptors_funcs

def _add_3d_coordinates(mol):
  mol = Chem.AddHs(mol)
  AllChem.EmbedMolecule(mol, AllChem.ETKDG())
  mol = Chem.RemoveHs(mol)
  return mol

_descriptors = list(descriptors)

length = len(smiles)
values = [np_none(length) for _ in descriptors]
descriptors_3d = ['PMI1', 'PMI2', 'PMI3', 'NPR1', 'NPR2', 'RadiusOfGyration',
                  'InertialShapeFactor', 'Eccentricity', 'Asphericity', 'SpherocityIndex']
descriptors_funcs = _get_descriptors_funcs(_descriptors)
for n in range(0, length):
  mol = Chem.MolFromMolBlock(smiles[n], sanitize = True) if ("M  END" in smiles[n]) else Chem.MolFromSmiles(smiles[n], sanitize = True)

  if mol is None:
    continue
  try:
    for d in _descriptors:
      if d in descriptors_3d:
        mol = _add_3d_coordinates(mol)
        break
    for d in range(0, len(_descriptors)):
      value = descriptors_funcs[_descriptors[d]](mol)
      if not np.isnan(value) and not np.isinf(value):
        values[d][n] = value
  except:
    continue
result = {_descriptors[d]: values[d] for d in range(len(_descriptors))}

out = pd.DataFrame(result)