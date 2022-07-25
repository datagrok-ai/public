#name: Embed
#language: python
#input: string molecule
#output: string sdf

from rdkit.Chem import AllChem
from rdkit import Chem
mol = AllChem.MolFromMolBlock(molecule) if ("M  END" in molecule) else AllChem.MolFromSmiles(molecule)

AllChem.EmbedMolecule(mol, AllChem.ETKDG())
#AllChem.UFFOptimizeMolecule(mol)
#mol = Chem.RemoveHs(mol)
sdf = Chem.MolToMolBlock(mol)
