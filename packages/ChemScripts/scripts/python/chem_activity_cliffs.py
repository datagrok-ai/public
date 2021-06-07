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
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.manifold import TSNE

similarityValue = 50
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
npfingerprints = []

for i in range(0, length):
    mol = Chem.MolFromSmiles(smiles[i], sanitize=True)
    if mol is None or mol.GetNumAtoms() == 0:
        continue
    mols[i] = mol
    accepted.append(i)

for i in accepted:
    fingerprints.append(Chem.RDKFingerprint(mols[i]))
    arr = np.zeros((0,))
    fp = AllChem.GetMorganFingerprintAsBitVect(mols[i], 2)
    DataStructs.ConvertToNumpyArray(fp, arr)
    npfingerprints.append(arr)

output_pairs = pd.DataFrame(columns=['n1', 'n2', "sim", "difference", "sali"])
maxSALI = 0
realLength = len(accepted)

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
        output_pairs.loc[i] = [np.int16(i), np.int16(j + i + 1), sim, activityDif, saliVal]

optSimilarityLimit = initialSimilarityLimit

output_pairs.loc[output_pairs['sim'] > optSimilarityLimit]

neighboursCount = [0] * length
similarityCount = [0] * length
separatedGraphs = [[]]
groupNumbers = [0] * length

for i in range(0, len(output_pairs.index)):
    neighboursCount[output_pairs.loc[i, "n1"].astype(int)] += 1
    neighboursCount[output_pairs.loc[i, "n2"].astype(int)] += 1
    similarityCount[output_pairs.loc[i, "n1"].astype(int)] += output_pairs.loc[i, "sim"]
    similarityCount[output_pairs.loc[i, "n2"].astype(int)] += output_pairs.loc[i, "sim"]

notAssigned = list(range(0, length))


def extract_sub_graph(node, used, pairs):
    if node in used:
        return list()

    neighbours = []

    for i in range(0, len(pairs.index)):
        if node == pairs.loc[i, "n1"].astype(int):
            neighbours.append(pairs.loc[i, "n2"].astype(int))
        elif node == pairs.loc[i, "n2"].astype(int):
            neighbours.append(pairs.loc[i, "n1"].astype(int))

    graph_nodes = [node]
    used.append(node)

    for i in range(0, len(neighbours)):
        new_nodes = extract_sub_graph(neighbours[i], used, pairs)
        graph_nodes += new_nodes
        used += new_nodes

    return graph_nodes


while len(notAssigned) > 0:
    graph = extract_sub_graph(notAssigned[0], list(), output_pairs)
    if len(graph) == 1:
        separatedGraphs[0] += graph
    else:
        separatedGraphs.append(graph)

    notAssigned = list(set(notAssigned) - set(graph))

for i in range(0, len(separatedGraphs)):
    for j in separatedGraphs[i]:
        groupNumbers[j] = i


def tanimoto_dist(a, b):
    dotprod = np.dot(a, b)
    tc = dotprod / (np.sum(a) + np.sum(b) - dotprod)
    return 1.0 - tc


tsne = TSNE(n_components=components, metric=tanimoto_dist)
tsne_X = tsne.fit_transform(npfingerprints)
x = tsne_X.T[0]
y = tsne_X.T[1]

output_coords = pd.DataFrame({'graph': groupNumbers, "x_coord": x, "y_coord": y})
