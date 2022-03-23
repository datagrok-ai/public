#name: Two Component Reaction
#description: Two component reaction
#help-url: https://datagrok.ai/help/domains/chem/functions/reactions
#language: python
#sample: chem/reactants.csv
#tags: demo, chem, rdkit
#input: dataframe data1 [First data table]
#input: column reactants1 {type:categorical; semType: Molecule} [Reactants molecules first set, in SMILES format]
#input: dataframe data2 [Second data table]
#input: column reactants2 {type:categorical; semType: Molecule} [Reactants molecules second set, in SMILES format]
#input: string reaction = [C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3] [Reaction, in SMARTS format]
#input: bool matrixExpansion = false [If checked, reactants will be combined from two sets, if is not - reactants will be combined sequentially]
#input: bool randomize = false [Randomize reactants]
#input: int seed = -1 [Random seed, set -1 to disable seed]
#input: int maxRandomReactions = 100 [Maximum number of reactions to be calculated]
#output: dataframe result {semType: Molecule} [Reactions results, in SMILES format]

import random
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

reactants1 = data1[reactants1]
reactants2 = data2[reactants2]
reaction = AllChem.ReactionFromSmarts(reaction)

def run_reaction(r1, r2):
    products = reaction.RunReactants((r1, r2))
    unique = {}
    for product in products:
        unique[Chem.MolToSmiles(product[0])] = products[0]
    return sorted(unique.keys())

column_names = ['reactant1', 'reactant2']
length1 = len(reactants1)
length2 = len(reactants2)
result_length = 0
max_rr_length = 0

if matrixExpansion:
    if randomize:
        if seed != -1:
            random.seed(seed)
        allocation_length = maxRandomReactions
    else:
        allocation_length = length1 * length2
else:
    allocation_length = min([length1, length2])

reactants1_ = np.full(allocation_length, None, dtype=object)
reactants2_ = np.full(allocation_length, None, dtype=object)
reaction_results = np.full(allocation_length, None, dtype=object)

def push_result(rr, reactant1_smiles, reactant2_smiles):
    global max_rr_length, result_length
    if len(rr) > max_rr_length:
        max_rr_length = len(rr)
    reaction_results[result_length] = rr
    reactants1_[result_length] = reactant1_smiles
    reactants2_[result_length] = reactant2_smiles
    result_length += 1

if matrixExpansion:
    if randomize:
        for n in range(0, length1 * length2):
            reactant1_smiles = reactants1[random.randint(0, length1 - 1)]
            reactant2_smiles = reactants2[random.randint(0, length2 - 1)]
            reactant1 = Chem.MolFromSmiles(reactant1_smiles)
            reactant2 = Chem.MolFromSmiles(reactant2_smiles)
            if reactant1 is None or reactant2 is None:
                continue
            rr = run_reaction(reactant1, reactant2)
            if len(rr) == 0:
                continue
            push_result(rr, reactant1_smiles, reactant2_smiles)
            if result_length >= maxRandomReactions:
                break
    else:
        for n in range(0, length1):
            reactant1_smiles = reactants1[n]
            reactant1 = Chem.MolFromSmiles(reactant1_smiles)
            for m in range(0, length2):
                reactant2_smiles = reactants2[m]
                reactant2 = Chem.MolFromSmiles(reactant2_smiles)
                if reactant1 is None or reactant2 is None:
                    continue
                rr = run_reaction(reactant1, reactant2)
                if len(rr) == 0:
                    continue
                push_result(rr, reactant1_smiles, reactant2_smiles)
else:
    for n in range(0, allocation_length):
        reactant1_smiles = reactants1[n]
        reactant2_smiles = reactants2[n]
        reactant1 = Chem.MolFromSmiles(reactants1[n])
        reactant2 = Chem.MolFromSmiles(reactants2[n])
        if reactant1 is None or reactant2 is None:
            continue
        rr = run_reaction(reactant1, reactant2)
        if len(rr) == 0:
            continue
        push_result(rr, reactant1_smiles, reactant2_smiles)

result = pd.DataFrame(np.column_stack([reactants1_[:result_length], reactants2_[:result_length]]), columns=column_names)
reaction_results_ = np.full((max_rr_length, result_length), None, dtype=object)
for n in range(0, result_length):
    rr = reaction_results[n]
    if rr is not None:
        for m in range(0, len(rr)):
            reaction_results_[m][n] = rr[m]
for m in range(0, max_rr_length):
    result['ps' + str(m + 1)] = reaction_results_[m]
