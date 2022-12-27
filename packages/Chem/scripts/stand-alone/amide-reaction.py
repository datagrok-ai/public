#name: Amide reaction
#language: python
#sample: chem/amines.csv, chem/carb_acids.csv
#tags: demo, chem, rdkit
#input: dataframe amines [Input data table]
#input: column amine_molecules {semType: Molecule}
#input: dataframe acids [Input data table]
#input: column acid_molecules {semType: Molecule}
#output: dataframe productsDf

from rdkit import Chem
from rdkit.Chem import rdChemReactions

rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')

def getMol(m):
    return Chem.MolFromSmiles(m, sanitize = True) if m is not None and "M  END" not in m else Chem.MolFromMolBlock(m, sanitize = True)

def getMolsArray(col):
    molsArray = []
    for molecule in col:
        molsArray.append(getMol(molecule))
    return molsArray
  
def getSmilesArray(col):
    smilesArray = []
    for molecule in col:
        smilesArray.append(Chem.MolToSmiles(molecule))
    return smilesArray


#convert smiles or molfiles to rdkit objects
amineMols = getMolsArray(amines[amine_molecules])
acidMols = getMolsArray(acids[acid_molecules])

amineCol = []
acidCol = []
productsArray = []

#create Cartesian product of reactants
for amIdx, am in enumerate(amineMols):
    for acIdx, ac in enumerate(acidMols):
        products = rxn.RunReactants((ac, am))
        productsArray.append(Chem.MolToSmiles(products[0][0])) if len(products) and len(products[0]) else productsArray.append('')
        amineCol.append(amines[amine_molecules][amIdx])
        acidCol.append(acids[acid_molecules][acIdx])


# Convert to Pandas DataFrame
aminesDf = pd.DataFrame(amineCol, columns=['amine'])
acidsDf = pd.DataFrame(acidCol, columns=['acid'])
prodDf = pd.DataFrame(productsArray, columns=['product'])
productsDf = pd.concat([aminesDf, acidsDf, prodDf], axis=1)