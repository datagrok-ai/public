#name: Similarity Maps Using Fingerprints
#description: Similarity Maps Using Fingerprints, RDKit based
#help-url: https://datagrok.ai/help/domains/chem/functions/sim-maps
#language: python
#tags: demo, chem, rdkit
#input: string mol = COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21 {semType: Molecule} [Molecule, in SMILES format]
#input: string refmol = CCCN(CCCCN1CCN(c2ccccc2OC)CC1)Cc1ccc2ccccc2c1 {semType: Molecule} [Reference molecule, in SMILES format]
#input: int radius = 2 [Fingerprint function radius]
#output: double maxweight [The maximum weight that was found when creating the map]
#output: graphics simMap [Similarity Map]

from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps

mol = Chem.MolFromSmiles(mol)
refmol = Chem.MolFromSmiles(refmol)

from rdkit import DataStructs
simMap, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, lambda m,
    idx: SimilarityMaps.GetMorganFingerprint(m, atomId=idx, radius=radius, fpType='count'),
    metric=DataStructs.TanimotoSimilarity)
