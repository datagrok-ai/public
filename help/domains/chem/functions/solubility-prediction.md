<!-- TITLE: Solubility prediction -->
<!-- SUBTITLE: -->

# Solubility prediction

Solubility is one of basic physical chemistry properties important for understanding how molecules interact with
solvents. Following script allows to predict solubility by molecular
descriptors. `#{x.18b704d0-0b50-11e9-b846-1fa94a4da5d1."Predict Solubility"}` model was trained
using `#{x.Demo:SolubilityTrain."Solubility Train"}` dataset on H2O modelling engine. Modelling method is "Generalized
Linear Modeling"

Molecular Descriptors used in model:

* **MolWt** - Molecular Weight
* **Ipc** - The information content of the coefficients of the characteristic polynomial of the adjacency matrix of a
  hydrogen-suppressed graph of a molecule
* **TPSA** - Total Polar Surface Area
* **LabuteASA** - Labute's Approximate Surface Area
* **NumHDonors** - Number of Hydrogen Donors
* **NumHAcceptors** - Number of Hydrogen Acceptors
* **MolLogP** - Wildman-Crippen LogP value
* **HeavyAtomCount** - Number of Heavy Atoms
* **NumRotatableBonds** - Number of Rotatable Bonds
* **RingCount** - Number of Rings
* **NumValenceElectrons** - Number of Valence Electrons

See also:

* [Cheminformatics](../cheminformatics.md)

References:

* [RDKit](https://www.rdkit.org)
* [MoleculeNet: A benchmark for molecular machine learning](https://arxiv.org/abs/1703.00564)
