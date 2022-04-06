#name: Gasteiger Partial Charges
#description: The Gasteiger partial charges visualization, RDKit based
#help-url: https://datagrok.ai/help/domains/chem/functions/gasteiger-charges
#language: python
#tags: demo, chem, rdkit, panel
#condition: true
#input: string mol = "COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21" {semType: Molecule} [Molecule, in SMILES format]
#input: int contours = 10
#output: graphics charges [The Gasteiger partial charges]

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps

mol = Chem.MolFromMolBlock(mol) if ("M  END" in mol) else Chem.MolFromSmiles(mol)
if mol is not None:
    AllChem.ComputeGasteigerCharges(mol)
    contribs = [float(mol.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(mol.GetNumAtoms())]
    charges = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, contourLines=contours)
