#name: FindRGroupsWithCore
#language: python
#input: string molecules
#input: dataframe df
#input: string core
#input: bool onlyMatchAtRGroups
#output: dataframe result

from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition
import numpy as np
import re

molecules = df[molecules].tolist()
length = len(molecules)
mols = []
r_group_map = dict()
ps = rdRGroupDecomposition.RGroupDecompositionParameters()
ps.onlyMatchAtRGroups = onlyMatchAtRGroups
core = Chem.MolFromMolBlock(core, sanitize = True) if ("M  END" in core) else Chem.MolFromSmarts(core)
if core is not None:
    for n in range(0, length):
        try:
          mol = Chem.MolFromMolBlock(molecules[n], sanitize = True) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n], sanitize = True)
        except:
          mol = None
        if mol is None:
          continue
        mols.append(mol)
    rgd,fails = rdRGroupDecomposition.RGroupDecompose([core], mols, asRows = False, options = ps)

    for key in rgd.keys():
        r_group_map[key] = []
        fail_counter = 0
        for j in range(0, length):
            if key == 'Core':
                r_group_map[key].append(Chem.MolToSmiles(rgd[key][0]))
            else:
            	if len(fails) > fail_counter and fails[fail_counter] == j:
                	fail_counter += 1
                	r_group_map[key].append('')
            	else:
              		r_group_map[key].append(Chem.MolToSmiles(rgd[key][j - fail_counter]))
else:
  raise Exception('Core is empty')                        
result = pd.DataFrame.from_dict(r_group_map)