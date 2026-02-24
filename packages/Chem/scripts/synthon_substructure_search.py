#name: SynthonSubstructureSearch
#description: Substructure search in synthon chemical space using RDKit SynthonSpaceSearch
#language: python
#environment: channels: [conda-forge], dependencies: [python=3.10, rdkit=2025.09.3]
#meta.cache: all
#meta.cache.invalidateOn: 0 0 * * *
#input: string molecule {semType: Molecule} [Query molecule in SMILES or Molblock format]
#input: file synthonLibrary
#input: string libraryName [Original library file name for caching]
#input: int maxHits = 100 [Maximum number of hit molecules to return]
#output: dataframe result

import os
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdSynthonSpaceSearch, rdFingerprintGenerator

if 'M  END' in molecule:
  mol = Chem.MolFromMolBlock(molecule, sanitize=True)
else:
  mol = Chem.MolFromSmiles(molecule, sanitize=True)

if mol is None:
  result = pd.DataFrame({'smiles': pd.Series(dtype='str'), 'name': pd.Series(dtype='str'), 'reaction_id': pd.Series(dtype='str'), 'component_count': pd.Series(dtype='int')})
else:
  synthons_dir = 'synthons'
  base_name = os.path.splitext(libraryName)[0]
  cached_db = os.path.join(synthons_dir, base_name + '.db')

  print(f'*** Checking cached DB: {cached_db}, exists: {os.path.isfile(cached_db)}')
  time.sleep(1)

  if not os.path.isfile(cached_db):
    print('*** DB not found, creating...')
    time.sleep(1)
    os.makedirs(synthons_dir, exist_ok=True)
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
    ss = rdSynthonSpaceSearch.SynthonSpace()
    ss.ReadTextFile(synthonLibrary)
    ss.BuildSynthonFingerprints(fpgen)
    ss.WriteDBFile(cached_db)
    print(f'*** DB created, exists now: {os.path.isfile(cached_db)}')
    time.sleep(1)
  else:
    print('*** Using cached DB')
    time.sleep(1)

  synthonspace = rdSynthonSpaceSearch.SynthonSpace()
  synthonspace.ReadDBFile(cached_db)
  print('*** DB loaded successfully')
  time.sleep(1)

  params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
  params.maxHits = maxHits

  search_results = synthonspace.SubstructureSearch(mol, params=params)
  hits = search_results.GetHitMolecules()

  smiles_list = []
  names_list = []
  reaction_id_list = []
  component_count_list = []
  for h in hits:
    smiles_list.append(Chem.MolToSmiles(h))
    try:
      name = h.GetProp('_Name')
      names_list.append(name)
      parts = name.split(';')
      reaction_id_list.append(parts[-1] if len(parts) > 1 else '')
      component_count_list.append(len(parts) - 1 if len(parts) > 1 else 0)
    except:
      names_list.append('')
      reaction_id_list.append('')
      component_count_list.append(0)

  result = pd.DataFrame({'smiles': smiles_list, 'name': names_list, 'reaction_id': reaction_id_list, 'component_count': component_count_list})
