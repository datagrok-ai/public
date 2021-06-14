#name: ActivityCliffsPy
#language: python
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#input: column activities
#input: double similarityValue
#output: dataframe output_coords
#output: dataframe output_pairs

import numpy as np
import math
from rdkit import Chem
import random
from rdkit import DataStructs

automaticSimilarityLimit = False
components = 3

AVERAGE_NEIGHBOR_COUNT = 6
MIN_SIMILARITY = 80
DEFAULT_SIMILARITY = 95
VIEW_CYCLE_COUNT = 5000

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
saliCount = [0] * length

for i in range(0, len(output_pairs.index)):
    neighboursCount[output_pairs.loc[i, "n1"].astype(int)] += 1
    neighboursCount[output_pairs.loc[i, "n2"].astype(int)] += 1
    similarityCount[output_pairs.loc[i, "n1"].astype(int)] += output_pairs.loc[i, "sim"]
    similarityCount[output_pairs.loc[i, "n2"].astype(int)] += output_pairs.loc[i, "sim"]
    if output_pairs.loc[i, "sali"] != math.inf:
        if activities[output_pairs.loc[i, "n1"]] > activities[output_pairs.loc[i, "n2"]]:
            saliCount[output_pairs.loc[i, "n1"].astype(int)] += output_pairs.loc[i, "sali"]
        else:
            saliCount[output_pairs.loc[i, "n2"].astype(int)] += output_pairs.loc[i, "sali"]

mx = []
my = []
mDx = []
mDy = []


for i in range(0, length):
   mx.append(random.random())
   my.append(random.random())

sorted_index = np.argsort(mx)
rowCount = length
minDistance = 1/math.sqrt(rowCount)
for cycle in range(0, VIEW_CYCLE_COUNT):
    mDx = [0] * rowCount
    mDy = [0] * rowCount

    cycleState = cycle / VIEW_CYCLE_COUNT
    cycleMinDistance = cycleState * (2 - cycleState) * minDistance
    if cycleState < 0.5:
        attractionCycleFactor = 0.8 * (1 - cycleState) * 0.5
        repulsionCycleFactor = 0.5
    else:
        attractionCycleFactor = 0.8 * (1 - cycleState) * (1 - cycleState)
        repulsionCycleFactor = 1 - cycleState

    for i in range(0, len(output_pairs.index)):
        n1 = output_pairs.loc[i, "n1"].astype(int)
        n2 = output_pairs.loc[i, "n2"].astype(int)
        dx = mx[n2] - mx[n1]
        dy = my[n2] - my[n1]
        distance = math.sqrt(dx*dx + dy*dy)
        shift = distance - minDistance
        neigbourFactor1 = 0
        neigbourFactor2 = 0
        if shift > 0:
            dx *= shift / distance
            dy *= shift / distance
            if neighboursCount[n1] > 4:
                neigbourFactor1 = 4.0/neighboursCount[n1]
            else:
                neigbourFactor1 = 1.0
            if neighboursCount[n2] > 4:
                neigbourFactor2 = 4.0/neighboursCount[n2]
            else:
                neigbourFactor2 = 1

            mDx[n1] = mDx[n1] + dx * attractionCycleFactor * neigbourFactor1
            mDy[n1] = mDy[n1] + dy * attractionCycleFactor * neigbourFactor1
            mDx[n2] = mDx[n2] - dx * attractionCycleFactor * neigbourFactor2
            mDy[n2] = mDy[n2] - dy * attractionCycleFactor * neigbourFactor2

    for i in range(0, len(sorted_index)):
        for j in range(i+1, len(sorted_index)):
            n1 = sorted_index[i]
            n2 = sorted_index[j]
            dx = mx[n2] - mx[1]
            if abs(dx) >= cycleMinDistance:
                break
            dy = my[n2] - my[n1]
            if abs(dy) < cycleMinDistance:
                distance = math.sqrt(dx * dx + dy * dy)
                if distance < cycleMinDistance:
                    if distance == 0:
                        angle = 2 * math.pi * random.random()
                        dx = math.sin(angle) * cycleMinDistance
                        dy = math.cos(angle) * cycleMinDistance
                    else:
                        shift = cycleMinDistance - distance
                        dx *= shift / distance
                        dy *= shift / distance

                    mDx[n1] = mDx[n1] - dx * repulsionCycleFactor
                    mDy[n1] = mDy[n1] - dy * repulsionCycleFactor
                    mDx[n2] = mDx[n2] + dx * repulsionCycleFactor
                    mDy[n2] = mDy[n2] + dy * repulsionCycleFactor

    for i in range(0, rowCount):
        mx[i] = min(max(mx[i] + mDx[i], 0), 1)
        my[i] = min(max(my[i] + mDy[i], 0), 1)

vec_lengths = [0]*rowCount
for i in range(0, rowCount):
    vec_lengths[i] = math.sqrt((mx[i]-0.5)*(mx[i]-0.5) + (my[i]-0.5)*(my[i]-0.5))

sorted_index = np.argsort(vec_lengths)

for i in range(0, rowCount):
    n = sorted_index[i]
    x = mx[n] - 0.5
    y = my[n] - 0.5
    a = 0
    if y != 0:
        a = math.atan(x / y)
        if y < 0:
            if x < 0:
                a -= math.pi
            else:
                a += math.pi
    else:
        if x > 0:
            a = math.pi / 2
        else:
            a = -math.pi / 2

    mx[n] = mx[n]*2 - 1
    my[n] = my[n]*2 - 1

output_coords = pd.DataFrame({"n": range(1, rowCount + 1), "x_coord": mx, "y_coord": my, "sali": saliCount})
