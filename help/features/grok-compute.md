<!-- TITLE: GrokCompute -->
<!-- SUBTITLE: -->

# GrokCompute

GrokCompute is high performance computation server. It is distributed as part of [Compute VM](compute-vm.md) 
that provided as [docker container](https://www.docker.com/).

## REST API

GrokCompute exposes API for the following features:

* Cheminformatics (RDKit-based)
    - Substructure search
    - Descriptors
    - Parse SDF
    - Find MCS
    - Get R-Groups
    - Smiles to 3d coordinates
    - Smiles to fingerprints
    - Smiles to InChI
    - Smiles to InChI key
    - InChI to InChI key
    - InChI to smiles
    - Smiles to Canonical
    - Draw molecule
    - Draw reaction
* Utilities
    - Jupyter Notebooks converter (HTML, PDF etc.)
    - PDF to CSV converter
* Modeling (concept only)
    - Train
    - Apply


See also:
* [Compute VM](compute-vm.md)
* [Cheminformatics](../domains/chem/cheminformatics.md)
