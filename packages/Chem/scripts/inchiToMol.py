#name: inchiToMol
#language: python
#meta.role: converter
#meta.inputRegexp: (InChI\=.+)
#connection: Chem
#input: string id
#output: string smiles { semType: Molecule }

from rdkit import Chem    
smiles = Chem.MolToSmiles(Chem.MolFromInchi(id))
