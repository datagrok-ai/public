#name: SmilesTo3DCoordinates
#language: python
#input: string molecule
#output: string sdf

from rdkit.Chem import AllChem
from rdkit import Chem

mol = Chem.MolFromSmiles(molecule)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
sdf = Chem.MolToMolBlock(mol)
