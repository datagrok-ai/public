#name: smilesTo3DCoordinates
#language: python
#input: string smiles
#output: string sdf

from rdkit.Chem import AllChem
from rdkit import Chem

mol = Chem.MolFromSmiles(smiles)
mol.SetProp("_Name", smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
mol = Chem.RemoveHs(mol)
sdf = Chem.MolToMolBlock(mol)
