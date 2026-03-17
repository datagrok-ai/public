#name: SynthonSearch
#description: Search in synthon chemical space using RDKit SynthonSpaceSearch
#language: python
#environment: channels: [conda-forge], dependencies: [python=3.10, rdkit=2025.09.3]
#meta.cache: all
#meta.cache.invalidateOn: 0 0 * * *
#input: string molecule {semType: Molecule} [Query molecule in SMILES or Molblock format]
#input: file synthonLibrary
#input: string libraryName [Original library file name for caching]
#input: int maxHits = 100 [Maximum number of hit molecules to return]
#input: string searchType = "substructure" {choices: ["substructure", "similarity", "exact"]} [Search type]
#input: double similarityCutoff = 0.5 [Minimum Tanimoto similarity threshold (used for similarity/exact)]
#input: bool includeSynthons = false [Include synthon structures and IDs]
#output: dataframe result

import os
import csv
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdSynthonSpaceSearch, rdFingerprintGenerator

if 'M  END' in molecule:
  mol = Chem.MolFromMolBlock(molecule, sanitize=True)
else:
  if '#' in molecule:
    mol = Chem.MolFromSmarts(molecule)
    if mol is None:
      mol = Chem.MolFromSmiles(molecule, sanitize=True)
  else:
    mol = Chem.MolFromSmiles(molecule, sanitize=True)

is_exact = searchType == 'exact'
is_similarity = searchType == 'similarity' or is_exact
include_similarity = is_similarity and not is_exact
cutoff = 1.0 if is_exact else similarityCutoff

if mol is None:
  result = pd.DataFrame()
else:
  synthons_dir = 'synthons'
  base_name = os.path.splitext(libraryName)[0]
  cached_db = os.path.join(synthons_dir, base_name + '.db')

  print(f'*** Checking cached DB: {cached_db}, exists: {os.path.isfile(cached_db)}')

  if not os.path.isfile(cached_db):
    print('*** DB not found, creating...')
    os.makedirs(synthons_dir, exist_ok=True)
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
    ss = rdSynthonSpaceSearch.SynthonSpace()
    ss.ReadTextFile(synthonLibrary)
    ss.BuildSynthonFingerprints(fpgen)
    ss.WriteDBFile(cached_db)
    print(f'*** DB created, exists now: {os.path.isfile(cached_db)}')
  else:
    print('*** Using cached DB')

  synthonspace = rdSynthonSpaceSearch.SynthonSpace()
  synthonspace.ReadDBFile(cached_db)
  print('*** DB loaded successfully')

  params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
  params.maxHits = maxHits

  t_search_start = time.time()
  if is_similarity:
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
    params.similarityCutoff = cutoff
    search_results = synthonspace.FingerprintSearch(mol, fpgen, params)
  else:
    search_results = synthonspace.SubstructureSearch(mol, params=params)

  hits = search_results.GetHitMolecules()
  t_search_end = time.time()
  print(f'*** Search completed in {t_search_end - t_search_start:.3f}s, {len(hits)} hits')

  t_post_start = time.time()
  smiles_list = []
  names_list = []
  similarity_list = []
  for h in hits:
    smiles_list.append(Chem.MolToSmiles(h))
    try:
      names_list.append(h.GetProp('_Name'))
    except:
      names_list.append('')
    if include_similarity:
      try:
        similarity_list.append(float(h.GetProp('Similarity')))
      except:
        similarity_list.append(0.0)

  columns = {'smiles': smiles_list}
  if include_similarity:
    columns['similarity'] = similarity_list

  if not includeSynthons:
    result = pd.DataFrame(columns)
  else:
    # Build synthon lookup from the library CSV
    synthon_lookup = {}
    with open(synthonLibrary, 'r', newline='') as f:
      reader = csv.DictReader(f)
      # Find the molecule column (SMILES) and ID column
      headers = reader.fieldnames or []
      smiles_col = None
      id_col = None
      for h_name in headers:
        lower = h_name.lower()
        if lower == 'smiles' or lower == 'molecule':
          smiles_col = h_name
        if lower == 'syntnon #' or lower == 'synton_id':
          id_col = h_name
      if smiles_col and id_col:
        for row in reader:
          synthon_lookup[row[id_col]] = row[smiles_col]

    # Parse names to extract synthon IDs and reaction IDs
    max_synthons = 0
    parsed = []
    for name in names_list:
      parts = name.split(';')
      if len(parts) > 1:
        synthon_ids = [s.strip() for s in parts[:-1]]
        reaction_id = parts[-1]
        parsed.append((reaction_id, synthon_ids))
        if len(synthon_ids) > max_synthons:
          max_synthons = len(synthon_ids)
      else:
        parsed.append(('', []))

    columns['reaction_id'] = [p[0] for p in parsed]
    for i in range(max_synthons):
      col_name = f'synthon_{i + 1}'
      col_id_name = f'synthon_{i + 1}_id'
      columns[col_name] = [
        synthon_lookup.get(p[1][i], '') if i < len(p[1]) else '' for p in parsed
      ]
      columns[col_id_name] = [
        p[1][i] if i < len(p[1]) else '' for p in parsed
      ]

    result = pd.DataFrame(columns)

  t_post_end = time.time()
  print(f'*** Post-processing completed in {t_post_end - t_post_start:.3f}s')
  print(f'*** Total: search {t_search_end - t_search_start:.3f}s + post-processing {t_post_end - t_post_start:.3f}s = {t_post_end - t_search_start:.3f}s')
