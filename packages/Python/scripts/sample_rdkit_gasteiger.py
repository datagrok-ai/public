#name: Sample Gasteiger charges 
#description: RDKit-based script.
#language: python
#input: string mol = "COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21" {semType: Molecule} [Molecule, in SMILES format]
#input: int contours = 10
#output: graphics charges [The Gasteiger partial charges]
#meta.queueName: python_docker

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps

global mol1
mol1 = Chem.MolFromMolBlock(mol) if ("M  END" in mol) else Chem.MolFromSmiles(mol)
if mol1 is not None:
	AllChem.ComputeGasteigerCharges(mol1)
	total = mol1.GetNumAtoms()
	contribs = [float(mol1.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(total)]
	charges = SimilarityMaps.GetSimilarityMapFromWeights(mol1, contribs, contourLines=contours)
	