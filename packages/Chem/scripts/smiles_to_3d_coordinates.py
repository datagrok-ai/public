#name: SmilesTo3DCoordinates
#language: python
#input: string smiles
#output: string sdf

from rdkit.Chem import AllChem
from rdkit import Chem

mol = Chem.MolFromSmiles(smiles)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
sdf = Chem.MolToMolBlock(mol)
