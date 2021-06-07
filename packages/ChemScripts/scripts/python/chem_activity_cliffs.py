#name: ActivityCliffsPy
#language: python
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#input: column activities
#output: dataframe output_coords
#output: dataframe output_pairs

import numpy as np
import math
from rdkit import Chem
import random
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.manifold import TSNE

similarityValue = 97
automaticSimilarityLimit = False
components = 3

AVERAGE_NEIGHBOR_COUNT = 6
MIN_SIMILARITY = 80
DEFAULT_SIMILARITY = 95
VIEW_CYCLE_COUNT = 20000

if automaticSimilarityLimit:
    initialSimilarityLimit = MIN_SIMILARITY
else:
    initialSimilarityLimit = similarityValue / 100

smiles = data[smiles]
activities = data[activities]
length = len(smiles)
mols = np.full(length, None, dtype=object)
accepted = []
fingerprints = []

for i in range(0, length):
    mol = Chem.MolFromSmiles(smiles[i], sanitize=True)
    if mol is None or mol.GetNumAtoms() == 0:
        continue
    mols[i] = mol
    accepted.append(i)

for i in accepted:
    fingerprints.append(Chem.RDKFingerprint(mols[i]))

pairs = []
maxSALI = 0
realLength = len(accepted)
pairs_count = 0
for i in range(0, realLength - 1):
    for j in range(0, realLength - i - 1):
        sim = DataStructs.TanimotoSimilarity(fingerprints[accepted[i]], fingerprints[accepted[j + i + 1]])
        if sim != 1:
            activityDif = abs(activities[accepted[i]] - activities[accepted[j + i + 1]])
            saliVal = activityDif / (1 - sim)
            if saliVal > maxSALI:
                maxSALI = saliVal
        else:
            saliVal = math.inf
            activityDif = 0
        pair = (i, j + i + 1, sim, activityDif, saliVal)
        pairs.append(pair)

optSimilarityLimit = initialSimilarityLimit

output_pairs = pd.DataFrame(pairs, columns=['n1', 'n2', "sim", "difference", "sali"])
output_pairs = output_pairs.loc[output_pairs['sim'] > optSimilarityLimit]
output_pairs = output_pairs.reset_index(drop=True)

neighboursCount = [0] * length
similarityCount = [0] * length
separatedGraphs = [[]]
groupNumbers = [0] * length

for i in range(0, len(output_pairs.index)):
    neighboursCount[output_pairs.loc[i, "n1"].astype(int)] += 1
    neighboursCount[output_pairs.loc[i, "n2"].astype(int)] += 1
    similarityCount[output_pairs.loc[i, "n1"].astype(int)] += output_pairs.loc[i, "sim"]
    similarityCount[output_pairs.loc[i, "n2"].astype(int)] += output_pairs.loc[i, "sim"]

mx = []
my = []

for i in range(0, length):
   mx.append(random.random())
   my.append(random.random())

output_coords = pd.DataFrame({"smiles": smiles, "activities":activities, "x_coord": mx, "y_coord": my})
