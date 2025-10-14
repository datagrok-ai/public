#name: Embed
#language: python
#input: string molecule
#output: string sdf

from rdkit.Chem import AllChem
from rdkit import Chem
mol = AllChem.MolFromMolBlock(molecule) if ("M  END" in molecule) else AllChem.MolFromSmiles(molecule)

sdf = Chem.MolToMolBlock(mol)
try:
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    sdf = Chem.MolToMolBlock(mol)
except Exception as e:
    pass
# mol = Chem.RemoveHs(mol)

