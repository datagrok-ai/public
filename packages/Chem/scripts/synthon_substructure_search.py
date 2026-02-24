#name: SynthonSubstructureSearch
#description: Substructure search in synthon chemical space using RDKit SynthonSpaceSearch
#language: python
#environment: channels: [conda-forge], dependencies: [python=3.9, rdkit=2024.9.3]
#meta.cache: all
#meta.cache.invalidateOn: 0 0 * * *
#input: string molecule {semType: Molecule} [Query molecule in SMILES or Molblock format]
#input: file synthonLibrary
#input: int maxHits = 100 [Maximum number of hit molecules to return]
#output: dataframe result

import os
import shutil
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdSynthonSpaceSearch, rdFingerprintGenerator

if 'M  END' in molecule:
  mol = Chem.MolFromMolBlock(molecule, sanitize=True)
else:
  mol = Chem.MolFromSmiles(molecule, sanitize=True)

if mol is None:
  result = pd.DataFrame({'smiles': pd.Series(dtype='str'), 'name': pd.Series(dtype='str')})
else:
  synthons_dir = 'synthons'
  file_name = os.path.basename(synthonLibrary)
  base_name = os.path.splitext(file_name)[0]
  cached_text = os.path.join(synthons_dir, file_name)
  cached_db = os.path.join(synthons_dir, base_name + '.db')

  if not os.path.isfile(cached_db):
    os.makedirs(synthons_dir, exist_ok=True)
    if not os.path.isfile(cached_text):
      shutil.copy2(synthonLibrary, cached_text)
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
    rdSynthonSpaceSearch.ConvertTextToDBFile(cached_text, cached_db, fpgen)

  synthonspace = rdSynthonSpaceSearch.SynthonSpace()
  synthonspace.ReadDBFile(cached_db)

  params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
  params.maxHits = maxHits

  search_results = synthonspace.SubstructureSearch(mol, params)
  hits = search_results.GetHitMolecules()

  smiles_list = []
  names_list = []
  for h in hits:
    smiles_list.append(Chem.MolToSmiles(h))
    try:
      names_list.append(h.GetProp('_Name'))
    except:
      names_list.append('')

  result = pd.DataFrame({'smiles': smiles_list, 'name': names_list})
