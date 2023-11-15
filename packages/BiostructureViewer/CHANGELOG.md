# BiostructureViewer changelog

## 1.1.0 (WIP)

### Features

* View for imported structures allows display controls
* Display input field to open structure on empty Biostructure viewer.
* Add pdbqt parser, convert to PDB, tests
* Pdbqt parser for models (poses) and target, tests
* Pdbqt import handler
* Context menu Download, Copy for Molecule3D grid cell
* Add Molecule3D ligands for Biostructure viewer

### Bug fixes

* Fix caching for PdbGridCellRenderer
* Fix MolstarViewer for Molecule ligands, tests
* Fix pdbGridCellRenderer, test

## 1.0.11 (2023-07-24)

*Dependency: datgarok-api >= 1.15.2*

### Features

* `Molecule3D`, `PDB_ID` semantic type
* [Biostructure](../../help/visualize/viewers/biostructure) viewer (Mol* based)
* [NGL](../../help/visualize/viewers/ngl) viewer
* [Molecule3D](https://public.datagrok.ai/apps/Tutorials/Demo/Bioinformatics/Proteins) cell renderer
* Biostructure and NGL viewers support [Molecule](../../help/develop/domains/chem) column for ligands
* File handlers for file [formats supported by Molstar or NGL](../../help/access/files/supported-formats.md)
