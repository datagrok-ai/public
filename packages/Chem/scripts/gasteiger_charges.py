#name: Chemistry | Gasteiger Partial Charges
#description: RDKit-based script.
#help-url: https://datagrok.ai/help/domains/chem/functions/gasteiger-charges
#language: python
#meta.role: panel
#meta.domain: chem
#meta.cache: all
#meta.cache.invalidateOn: 0 0 * * *
#condition: true
#input: string mol = "COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21" {semType: Molecule} [Molecule, in SMILES format]
#input: int contours = 10
#output: graphics charges [The Gasteiger partial charges]

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D
from IPython.display import Image, display

mol = Chem.MolFromMolBlock(mol) if ("M  END" in mol) else Chem.MolFromSmiles(mol)
if mol is not None:
    AllChem.ComputeGasteigerCharges(mol)
    contribs = [float(mol.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(mol.GetNumAtoms())]
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
    SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, draw2d=drawer, contourLines=contours)
    drawer.FinishDrawing()
    display(Image(data=drawer.GetDrawingText()))
