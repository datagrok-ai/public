#name: InchiToMolPython
#description: The InChI's molecule (inchiToMol's RDKit side: the converter Chem:InchiToMol calls this script, or Crux in Crux mode)
#language: python
#connection: Chem
#input: string id
#output: string smiles { semType: Molecule }

from rdkit import Chem    
smiles = Chem.MolToSmiles(Chem.MolFromInchi(id))
