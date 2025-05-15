# Supported formats

Datagrok supports over 50 file formats, including domain-specific like SDF, FASTA, and others. It also understands how data is structured or represented within those files (e.g., assay plates within an XLSX, SDTM conventions within a CSV, and notations like SMILES and HELM).

## Structured and analytical data

### Tabular and semi-structured data

| Format            | Type              | Description                                   | Required plugins |
|------------------|-------------------|-----------------------------------------------|------------------|
| `.csv`             | File format       | Comma-separated values                         |                  |
| `.tsv`             | File format       | Tab-separated values                           |                  |
| `.txt`             | File format       | Plain text                                     |                  |
| `.xlsx`            | File format       | Microsoft Excel                                |                  |
| `.d42`             | File format       | Datagrok [project]                             |                  |
| `.json`            | File format       | JavaScript Object Notation                     |                  |
| `.xml `            | File format       | Extensible Markup Language                     |                  |
| `.HTML`            | File format       | HyperText Markup Language                      |                  |

### Statistical and scientific formats

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.rds`, `.rda`       | File format       | R serialized data                              |                  |
| `.mat`             | File format       | MATLAB data format                             |                  |
| `.h5`              | File format       | Hierarchical Data Format                       |                  |
| `.nc`, `.netcdf`     | File format       | Network Common Data Form                       |                  |
| `.sas7bdat`        | File format       | SAS database file                              |                  |
| `.kxl`             | File format       | KeyCreator eXtensible Language                 |                  |

### Columnar and database formats

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.parquet`         | File format       | Columnar storage format                        | [Arrow]          |
| `.feather`         | File format       | Arrow-based binary format                      | [Arrow]          |
| `.sqlite`          | File format       | SQLite database                                | [SQLite]         |

### Archives and packaging

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.zip`             | File format       | ZIP archive (auto-extracted for supported types)|                 |
| `.tar`             | File format       | Tape archive                                   |                  |
| `.gz`, `.gzip`       | File format       | Gzip-compressed files                          |                  |


## Notebooks and documents

### Notebooks and computations

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.ipynb`           | File format       | Jupyter Notebook                               | [Notebooks]      |
| `.ivp`             | File format       | Interactive differential equations             | [Diff Studio]    |

### Docs and reports

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.pdf`             | File format       | Portable Document Format                       | [File Editors]   |
| `.docx`            | File format       | Microsoft Word Document                        | [File Editors]   |


## Domain-specific formats

### Cheminformatics

#### Molecular structures

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.sdf`, `.sd `       | File format       | Structure-data file                            | [Chem]           |
| `.mol`             | File format       | MDL Molfile                                    | [Chem]           |
| `.mol2`            | File format       | SYBYL molecule representation                  | [Chem]           |
| `.smi`             | File format       | Structure-data file                            | [Chem]           |

### Bioinformatics

#### Macromolecules and sequences

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.fasta`           | File format       | FASTA format (fa, fna, ffn, faa, frn, fa, fst) | [Bio]            |

#### Phylogenetic trees

| Format           | Type              | Description                                    | Required plugins     |
|------------------|-------------------|------------------------------------------------|----------------------|
| `.nwk`, `.newick`    | File format       | Newick format for phylogenetic trees           | [PhyloTree Viewer]   |

#### 3D structures and density maps

##### Macromolecular structure formats

| Format           | Type              | Description                                    | Required plugins        |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.pdb`             | File format       | Protein Data Bank                              | [Biostructure Viewer]    |
| `.pdbqt`           | File format       | Protein Data Bank, Partial Charge, Atom Type   | [Biostructure Viewer]    |
| `.cif`             | File format       | (preview) Crystallographic Information File    | [Biostructure Viewer]    |
| `.cifCore`         | File format       | Crystallographic Information Core File         | [Biostructure Viewer]    |
| `.mcif`            | File format       | (preview) Macromolecular Crystallographic File | [Biostructure Viewer]    |
| `.mmcif`           | File format       | Macromolecular Crystallographic File           | [Biostructure Viewer]    |
| `.mmtf`            | File format       | Macromolecular Transmission Format             | [Biostructure Viewer]    |
| `.ent`             | File format       | (preview) Brookhaven PDB Molecule              | [Biostructure Viewer]    |

##### Electron density and volumetric data

| Format           | Type              | Description                                    | Required plugins        |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.ccp4`            | File format       | Collaborative Computational Project Number 4  | [Biostructure Viewer]    |
| `.brix`            | File format       | (preview) Electron Density Map                 | [Biostructure Viewer]    |
| `.mrc`             | File format       | (preview) Electron Density Map                 | [Biostructure Viewer]    |
| `.xplor`           | File format       | (preview) Electron Density Map                 | [Biostructure Viewer]    |
| `.cub`, `.cube`      | File format       | (preview) Orbital/Density Values on a 3D Grid | [Biostructure Viewer]    |
| `.dsn6 `           | File format       | (preview) Ribbons FRODO map                    | [Biostructure Viewer]    |

##### Simulation and topology formats

| Format           | Type              | Description                                    | Required plugins        |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.gro`             | File format       | GROMACS structure                              | [Biostructure Viewer]    |
| `.top `            | File format       | GROMACS topology                               | [Biostructure Viewer]    |
| `.prmtop`          | File format       | AMBER Parameter/Topology File                  | [Biostructure Viewer]    |
| `.parm7`           | File format       | (preview) AMBER Parameter/Topology File        | [Biostructure Viewer]    |
| `.psf `            | File format       | (preview) Protein Structure File               | [Biostructure Viewer]    |
| `.cns`             | File format       | Crystallography and NMR System                 | [Biostructure Viewer]    |

##### 3D mesh formats

| Format           | Type              | Description                                    | Required plugins        |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.obj `            | File format       | Wavefront 3D object                            | [Biostructure Viewer]    |
| `.ply`             | File format       | Polygon file format                            | [Biostructure Viewer]    |
| `.xyz`             | File format       | (preview) Cartesian Coordinate Format          | [Biostructure Viewer]    |

### Spectroscopy and NMR

| Format           | Type              | Description                                    | Required plugins        |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.jdx `            | File format       | Spectroscopy                                   | [Chem Spectra Viewer]    |
| `.jdx`, `.dx`, `.nmrium`| File format       | NMR                                            | [NMRium]                 |

### Biosignals

| Format           | Type              | Description                                                     | Required plugins     |
|------------------|-------------------|------------------------------------------------------------------|----------------------|
| `.edf`             | File format       | European Data Format for biosignals                             |                      |

### Clinical and regulatory data

| Format           | Type              | Description                                                     | Required plugins     |
|------------------|-------------------|------------------------------------------------------------------|----------------------|
| `.csv`             | File format       | SDTM-formatted datasets (auto-detected from folder structure)   | [Clinical Case]      |
| `.xlsx`            | Data structure       | Assay plate data (auto-detected in file)                        | [Curves]             |

### Geospatial formats

| Format           | Type              | Description                                    | Required plugins |
|------------------|-------------------|------------------------------------------------|------------------|
| `.kml`, `.kmz`       | File format       | Keyhole Markup Language                        | [GIS]            |
| `.geojson`         | File format       | GeoJSON (geospatial data)                      | [GIS]            |
| `.topojson`        | File format       | TopoJSON (topological geospatial data)         | [GIS]            |

## Molecular and macromolecular notations (used in strings or columns)

| Format           | Type      | Description                                    | Required plugins |
|------------------|-----------|------------------------------------------------|------------------|
| `SMILES`| Notation  | Linear representation of molecular structure   | [Chem] |
| `SMARTS` | Notation  | Linear representation of molecular substructure patterns | [Chem]  |
| `InChI`, `InChIKey`| Notation  | IUPAC-based textual identifiers for chemical compounds | [Chem]  |
| `HELM` | Notation  | Biopolymer notation for peptides and oligonucleotides | [Bio] |
| `SMIRKS` | Notation  | Linear representation of chemical reaction transforms   | [Chem] |
| `BILN` | Notation  | Linear notation for complex biopolymers   | [Bio]  |

[Project]: <../../datagrok/concepts/project/project.md>
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