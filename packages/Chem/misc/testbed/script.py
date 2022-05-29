from rdkit import Chem, RDConfig
from rdkit.Chem import rdSubstructLibrary
from time import sleep, perf_counter as pc
import csv

t0 = pc()
library = rdSubstructLibrary.SubstructLibrary()
molholder = rdSubstructLibrary.CachedSmilesMolHolder()
patternHolder = rdSubstructLibrary.PatternHolder()
for i, row in enumerate(csv.reader(open('chembl_100k.csv', 'r'))):
  if 30000 <= i and i <= 40000:
    idx = molholder.AddSmiles(row[0])
    idx2 = patternHolder.AddFingerprint(
      patternHolder.MakeFingerprint(Chem.MolFromSmiles(row[0])))
    assert idx == idx2
library = rdSubstructLibrary.SubstructLibrary(molholder, patternHolder)
print("Construct: " + str(round(pc() - t0, 2)))

searchFor = ['c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O']
for smiles in searchFor:
  t0 = pc()
  core = Chem.MolFromSmarts('CCCCOC')
  indices = library.GetMatches(core, useChirality = False) 
  print("Search for " + smiles + ": " + str(round(pc() - t0, 2)))