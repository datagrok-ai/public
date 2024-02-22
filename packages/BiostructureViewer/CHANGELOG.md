# BiostructureViewer changelog

## 1.0.25 (2024-01-29)

## Features

* MolstarViewer optimize buildView postponed, destroyView detach only

### Bug fixes

* Fix MolstarViewer tests
* Fix BiotrackViewer for arbitrary data frame error
* Fix MolstarViewer test error Unsupported file extension bcif
* Fix MolstarViewer for pdb_data test

## 1.0.24 (2024-01-24)

### Bug fixes

* Fix error on right click on row header,test

## 1.0.23 (2024-01-19)

### Features

* View for imported structures allows display controls
* Display input field to open structure on empty Biostructure viewer.
* Add pdbqt parser, convert to PDB, tests
* Pdbqt parser for models (poses) and target, tests
* Pdbqt import handler
* Context menu Download, Copy for Molecule3D grid cell
* Add Molecule3D ligands for Biostructure viewer
* Add PdbHelper converters parsePdbqt, molToPdb, pdbqtToMol with NGL
* Add NgViewer support formats other than PDB
* Add columns descriptions to Demo Docking Conformations
* Add layout to Demo Docking Conformations
* Add AutoDockService
* Add Molecule3D detector for pdbqt units
* Split to Docking package

### Bug fixes

* Fix caching for PdbGridCellRenderer
* Fix MolstarViewer for Molecule ligands, tests
* Fix pdbGridCellRenderer, test
* Fix NglGlDocService to not miss tasks, refactor, handle errors, restore by timeout count limit
* Fix MolstarViewer for dataFrame detach
* Fix MolstarViewer usages to dispose WebGL on view closing, Could not create WebGL program
* Fix NglViewer ensure create, waits first render
* Fix MolstarViewer ensure creat, wait canvas3dInit
* Fix AutoDock to support atoms of a ligand
* Add AutoDock test
* Fix MolstarViewer, NglViewer tests for IRenderer.awaitRendered
* Fix MolstarViewer for retaining ligand
* Add tests for pdbqt from autodock-gpu
* Add pdbqt tests and sample data
* Add func/BiostructureViewer.autoDockApp
* Fix Demo to check empty data
* Fix import pdbqt to open poses with receptor
* Fix pdbqt parser for an element and atom type
* Fix pdbqt parser to get autodock-gpu scores
* Fix PdbGridCellRenderer, NglGlDocService to not mix cells
* Fix MolstarViewer to handle ligandColumnName changed
* Add test MolstarViewer for pdb_data with cell Click
* Add test MolstarViewer for pdb_id, data provider
* Fix biostructure data provider to return JSON string
* Add test biostructure provider with binary data (bcif)

## 1.0.11 (2023-07-24)

*Dependency: datgarok-api >= 1.15.2*

### Features

* `Molecule3D`, `PDB_ID` semantic type
* [Biostructure](../../help/visualize/viewers/biostructure) viewer (Mol* based)
* [NGL](../../help/visualize/viewers/ngl) viewer
* [Molecule3D](https://public.datagrok.ai/apps/Tutorials/Demo/Bioinformatics/Proteins) cell renderer
* Biostructure and NGL viewers support [Molecule](../../help/develop/domains/chem) column for ligands
* File handlers for file [formats supported by Molstar or NGL](../../help/access/files/supported-formats.md)
