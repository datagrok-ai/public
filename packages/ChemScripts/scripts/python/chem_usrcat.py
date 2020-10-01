#name: USRCAT
#description: USRCAT - real-time ultrafast shape recognition with pharmacophoric constraints
#help-url: https://datagrok.ai/help/domains/chem/functions/usrcat
#language: python
#sample: chem/smiles_coordinates.csv
#tags: demo, chem, rdkit
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#output: graphics scores

import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdBase
from rdkit.Chem import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT
from rdkit.Chem import DataStructs

smiles = data[smiles]
mols = [Chem.MolFromSmiles(mol) for mol in smiles]
mols = [mol for mol in mols if mol is not None]
for mol in mols:
    AllChem.EmbedMolecule(mol,
                          useExpTorsionAnglePrefs=True,
                          useBasicKnowledge=True)
usrcats = [GetUSRCAT(mol) for mol in mols]
fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2) for mol in mols]
data = {"tanimoto": [], "usrscore": []}
for i in range(len(usrcats)):
    for j in range(i):
        tc = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        score = GetUSRScore(usrcats[i], usrcats[j])
        data["tanimoto"].append(tc)
        data["usrscore"].append(score)
df = pd.DataFrame(data)
sns.pairplot(df)
