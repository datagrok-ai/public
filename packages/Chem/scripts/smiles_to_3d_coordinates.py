#name: SmilesTo3DCoordinates
#language: python
#input: string molecule
#output: string sdf

from rdkit.Chem import AllChem
from rdkit import Chem

mol = Chem.MolFromMolBlock(molecule, sanitize = True) if ("M  END" in molecule) else Chem.MolFromSmiles(molecule, sanitize = True)

AllChem.EmbedMolecule(mol, AllChem.ETKDG())
sdf = Chem.MolToMolBlock(mol)
