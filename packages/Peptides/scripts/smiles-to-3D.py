#name: smiTo3D
#language: python
#input: string smiles
#output: string sdf

from rdkit.Chem import AllChem
from rdkit import Chem

mol = AllChem.MolFromSmiles(smiles)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
#AllChem.UFFOptimizeMolecule(mol)
#mol = Chem.RemoveHs(mol)
sdf = Chem.MolToMolBlock(mol)
