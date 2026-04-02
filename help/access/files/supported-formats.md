# Supported formats

Datagrok supports over 50 file formats, including domain-specific formats like SDF, FASTA, and others. Datagrok also understands how data is structured or represented within those files (e.g., assay plates within an XLSX, SDTM conventions within a CSV, and notations like SMILES and HELM).

## Structured and analytical data

### Tabular and semi-structured data

| Format  | Description                | Required plugins |
|---------|----------------------------|------------------|
| `.csv`  | Comma-separated values     | --               |
| `.tsv`  | Tab-separated values       | --               |
| `.txt`  | Plain text                 | --               |
| `.xlsx` | Microsoft Excel            | --               |
| `.d42`  | Datagrok [Dashboard]       | --               |
| `.json` | JavaScript Object Notation | --               |
| `.xml ` | Extensible Markup Language | --               |
| `.HTML` | HyperText Markup Language  | --               |

### Statistical and scientific formats

| Format           | Description                        | Required plugins |
|------------------|------------------------------------|------------------|
| `.rds`, `.rda`   | R serialized data                  | --               |
| `.mat`           | MATLAB array and matrix data       | --               |
| `.h5`            | Hierarchical data in HDF5          | --               |
| `.nc`, `.netcdf` | NetCDF scientific data file        | --               |
| `.sas7bdat`      | SAS database file                  | --               |
| `.kxl`           | KeyCreator engineering data format | --               |
| `.pzfx`          | Prism 5                            | Curves           |
| `.prism`         | Prism 10                           | Prism            |

### Columnar and database formats

| Format     | Description               | Required plugins |
|------------|---------------------------|------------------|
| `.parquet` | Columnar storage format   | [Arrow]          |
| `.feather` | Arrow-based binary format | [Arrow]          |
| `.sqlite`  | SQLite database           | [SQLite]         |

### Archives and packaging

| Format         | Description                                      | Required plugins |
|----------------|--------------------------------------------------|------------------|
| `.zip`         | ZIP archive (auto-extracted for supported types) | --               |
| `.tar`         | Tape archive                                     | --               |
| `.gz`, `.gzip` | Gzip-compressed files                            | --               |


## Notebooks and documents

### Notebooks and computations

| Format   | Description                        | Required plugins |
|----------|------------------------------------|------------------|
| `.ipynb` | Jupyter Notebook                   | [Notebooks]      |
| `.ivp`   | Interactive differential equations | [Diff Studio]    |

### Docs and reports

| Format  | Description              | Required plugins |
|---------|--------------------------|------------------|
| `.pdf`  | Portable Document Format | [File Editors]   |
| `.docx` | Microsoft Word document  | [File Editors]   |
| `.tex`  | LaTeX Source File        | [File Editors]   |


## Domain-specific formats

### Cheminformatics

#### Molecular structures

| Format              | Description                                  | Required plugins |
|---------------------|----------------------------------------------|------------------|
| `.sdf`, `.sd `      | Structure-data file                          | [Chem]           |
| `.mol`              | MDL Molfile                                  | [Chem]           |
| `.mol2`             | SYBYL molecule representation                | [Chem]           |
| `.smi`              | SMILES file format                           | [Chem]           |
| `SMILES`            | Linear representation of molecular structure | [Chem]           |
| `SMARTS`            | Substructure search patterns                 | [Chem]           |
| `InChI`, `InChIKey` | IUPAC identifiers                            | [Chem]           |
| `SMIRKS`            | Reaction transformation notation             | [Chem]           |

### Bioinformatics

#### Macromolecules and sequences

| Format   | Description                                              | Required plugins |
|----------|----------------------------------------------------------|------------------|
| `.fasta` | FASTA format (fa, fna, ffn, faa, frn, fa, fst)           | [Bio]            |
| `HELM`   | Biopolymer notation for peptides and oligonucleotides    | [Bio]            |
| `BILN`   | Linear notation for complex biopolymers                  | [Bio]            |

#### Phylogenetic trees

| Format            | Description                          | Required plugins   |
|-------------------|--------------------------------------|--------------------|
| `.nwk`, `.newick` | Newick format for phylogenetic trees | [PhyloTree Viewer] |

#### 3D structures and density maps

##### Macromolecular structure formats

| Format     | Description                                 | Required plugins      |
|------------|---------------------------------------------|-----------------------|
| `.pdb`     | Protein structure (coordinates only)        | [Biostructure Viewer] |
| `.pdbqt`   | PDB with charges and atom types             | [Biostructure Viewer] |
| `.cif`     | Crystallographic data                       | [Biostructure Viewer] |
| `.cifCore` | Core crystallographic data                  | [Biostructure Viewer] |
| `.mcif`    | Macromolecular crystallographic data        | [Biostructure Viewer] |
| `.mmcif`   | Macromolecular crystallographic data        | [Biostructure Viewer] |
| `.mmtf`    | Compressed macromolecular structure format  | [Biostructure Viewer] |
| `.ent`     | Legacy PDB structure file                   | [Biostructure Viewer] |

##### Electron density and volumetric data

| Format          | Description                       | Required plugins      |
|-----------------|-----------------------------------|-----------------------|
| `.ccp4`         | Electron density map              | [Biostructure Viewer] |
| `.brix`         | Electron density map              | [Biostructure Viewer] |
| `.mrc`          | Electron density volume           | [Biostructure Viewer] |
| `.xplor`        | Electron density map              | [Biostructure Viewer] |
| `.cub`, `.cube` | 3D orbital or charge density grid | [Biostructure Viewer] |
| `.dsn6 `        | Ribbons-format electron map       | [Biostructure Viewer] |

##### Simulation and topology formats

| Format    | Description                       | Required plugins      |
|-----------|-----------------------------------|-----------------------|
| `.gro`    | GROMACS structure file            | [Biostructure Viewer] |
| `.top `   | GROMACS topology definition       | [Biostructure Viewer] |
| `.prmtop` | AMBER force field topology        | [Biostructure Viewer] |
| `.parm7`  | AMBER topology file               | [Biostructure Viewer] |
| `.psf `   | CHARMM or NAMD structure topology | [Biostructure Viewer] |
| `.cns`    | NMR structure description for CNS | [Biostructure Viewer] |

##### 3D mesh formats

| Format  | Description                      | Required plugins      |
|---------|----------------------------------|-----------------------|
| `.obj ` | 3D mesh geometry                 | [Biostructure Viewer] |
| `.ply`  | Polygon mesh data                | [Biostructure Viewer] |
| `.xyz`  | Atomic coordinates in XYZ format | [Biostructure Viewer] |


### Bioassays and experimental data

| Format  | Description                              | Required plugins |
|---------|------------------------------------------|------------------|
| `.xlsx` | Assay plate data (auto-detected in file) | [Curves]         |

### Spectroscopy and NMR

| Format    | Description                              | Required plugins      |
|-----------|------------------------------------------|-----------------------|
| `.jdx `   | Spectroscopy data (JCAMP-DX)             | [Chem Spectra Viewer] |
| `.dx`     | NMR spectroscopy data (JCAMP-DX variant) | [NMRium]              |
| `.nmrium` | NMR experiment data for NMRIum viewer    | [NMRium]              |

### Biosignals

| Format | Description                                          | Required plugins |
|--------|------------------------------------------------------|------------------|
| `.edf` | Multichannel biosignal recordings (e.g., EEG, EMG)   | --               |

### Clinical and regulatory data

| Format | Description                                                   | Required plugins |
|--------|---------------------------------------------------------------|------------------|
| `.csv` | SDTM-formatted datasets (auto-detected from folder structure) | [Clinical Case]  |

### Geospatial formats

| Format      | Description                                 | Required plugins |
|-------------|---------------------------------------------|------------------|
| `.kml`      | Geospatial annotation file (Google Earth)   | [GIS]            |
| `.kmz`      | Compressed geospatial annotation file       | [GIS]            |
| `.geojson`  | Geospatial data encoded in JSON             | [GIS]            |
| `.topojson` | Compressed topological geospatial data      | [GIS]            |

## Images and diagrams

| Format          | Description                  | Required plugins |
|-----------------|------------------------------|------------------|
| `.jpg`, `.jpeg` | Photo, scan                  |                  |
| `.png`          | Screenshot, icon, UI graphic |                  |
| `.svg`          | Logo, diagram, UI asset      |                  |
| `.bmp`          | Raw pixel data               |                  |
| `.excalidraw`   | JSON-based vector drawings   | [Excalidraw]     |



[Dashboard]: <../../datagrok/concepts/project/dashboard.md>
[Chem]: <https://github.com/datagrok-ai/public/tree/master/packages/Chem#readme>
[Bio]: <https://github.com/datagrok-ai/public/tree/master/packages/Bio#readme>
[Biostructure Viewer]: <https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer#readme>
[GIS]: https://github.com/datagrok-ai/public/tree/master/packages/GIS#readme
[Arrow]: https://github.com/datagrok-ai/public/tree/master/packages/Arrow#readme
[PhyloTree Viewer]: https://github.com/datagrok-ai/public/tree/master/packages/PhyloTreeViewer#readme
[File Editors]: https://github.com/datagrok-ai/public/tree/master/packages/FileEditors#readme
[Diff Studio]: https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio#readme
[SQLite]: https://github.com/datagrok-ai/public/tree/master/packages/SQLite#readme
[NMRium]: https://github.com/datagrok-ai/public/tree/master/packages/nmrium#readme
[Chem Spectra Viewer]: https://github.com/datagrok-ai/chem-spectra-viewer
[Notebooks]: https://github.com/datagrok-ai/public/tree/master/packages/Notebooks#readme
[Curves]: https://github.com/datagrok-ai/public/tree/master/packages/Curves#readme
[Clinical Case]: https://github.com/datagrok-ai/public/tree/master/packages/ClinicalCase#readme
[Excalidraw]: https://github.com/datagrok-ai/public/tree/master/packages/Excalidraw#readme
