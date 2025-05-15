# Supported formats

Datagrok supports over 50 file formats, including domain-specific like SDF, FASTA, and others. It also understands how data is structured or represented within those files (e.g., assay plates within an XLSX, SDTM conventions within a CSV, and notations like SMILES and HELM).

## Structured and analytical data

### Tabular and semi-structured data

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|-----------------------------------------------|------------------|
| `.csv`             | File format       | Comma-separated values                         | --               |
| `.tsv`             | File format       | Tab-separated values                           | --               |
| `.txt`             | File format       | Plain text                                     | --               |
| `.xlsx`            | File format       | Microsoft Excel                                | --               |
| `.d42`             | File format       | Datagrok [project]                             | --               |
| `.json`            | File format       | JavaScript Object Notation                     | --               |
| `.xml `            | File format       | Extensible Markup Language                     | --               |
| `.HTML`            | File format       | HyperText Markup Language                      | --               |

### Statistical and scientific formats

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.rds`, `.rda`       | File format       | R serialized data                            | --                 |
| `.mat`             | File format       | MATLAB array and matrix data                             | --                  |
| `.h5`              | File format       | Hierarchical data in HDF5                      | --                 |
| `.nc`, `.netcdf`     | File format       | NetCDF scientific data file                      | --                 |
| `.sas7bdat`        | File format       | SAS database file                              | --                 |
| `.kxl`             | File format       | KeyCreator engineering data format                 | --                 |

### Columnar and database formats

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.parquet`         | File format       | Columnar storage format                        | [Arrow]          |
| `.feather`         | File format       | Arrow-based binary format                      | [Arrow]          |
| `.sqlite`          | File format       | SQLite database                                | [SQLite]         |

### Archives and packaging

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.zip`             | File format       | ZIP archive (auto-extracted for supported types)| --                |
| `.tar`             | File format       | Tape archive                                   | --                 |
| `.gz`, `.gzip`       | File format       | Gzip-compressed files                          | --                 |


## Notebooks and documents

### Notebooks and computations

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.ipynb`           | File format       | Jupyter Notebook                               | [Notebooks]      |
| `.ivp`             | File format       | Interactive differential equations             | [Diff Studio]    |

### Docs and reports

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.pdf`             | File format       | Portable Document Format                       | [File Editors]   |
| `.docx`            | File format       | Microsoft Word document                        | [File Editors]   |


## Domain-specific formats

### Cheminformatics

#### Molecular structures

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.sdf`, `.sd `       | File format       | Structure-data file                            | [Chem]           |
| `.mol`             | File format       | MDL Molfile                                    | [Chem]           |
| `.mol2`            | File format       | SYBYL molecule representation                  | [Chem]           |
| `.smi`             | File format       | SMILES file format                           | [Chem]           |
| `SMILES`| Notation  | Linear representation of molecular structure   | [Chem] |
| `SMARTS` | Notation  | Substructure search patterns | [Chem]  |
| `InChI`, `InChIKey`| Notation  | IUPAC identifiers | [Chem]  |
| `SMIRKS` | Notation  | Reaction transformation notation   | [Chem] |

### Bioinformatics

#### Macromolecules and sequences

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.fasta`           | File format       | FASTA format (fa, fna, ffn, faa, frn, fa, fst) | [Bio]            |
| `HELM` | Notation  | Biopolymer notation for peptides and oligonucleotides | [Bio] |
| `BILN` | Notation  | Linear notation for complex biopolymers   | [Bio]  |

#### Phylogenetic trees

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|----------------------|
| `.nwk`, `.newick`    | File format       | Newick format for phylogenetic trees           | [PhyloTree Viewer]   |

#### 3D structures and density maps

##### Macromolecular structure formats

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.pdb`             | File format       | Protein structure (coordinates only)           | [Biostructure Viewer]    |
| `.pdbqt`           | File format       | PDB with charges and atom types                | [Biostructure Viewer]    |
| `.cif`             | File format       | Crystallographic data                          | [Biostructure Viewer]    |
| `.cifCore`         | File format       | Core crystallographic data                     | [Biostructure Viewer]    |
| `.mcif`            | File format       | Macromolecular crystallographic data | [Biostructure Viewer]    |
| `.mmcif`           | File format       | Macromolecular crystallographic data          | [Biostructure Viewer]    |
| `.mmtf`            | File format       | Compressed macromolecular structure format    | [Biostructure Viewer]    |
| `.ent`             | File format       | Legacy PDB structure file              | [Biostructure Viewer]    |

##### Electron density and volumetric data

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.ccp4`            | File format       | Electron density map  | [Biostructure Viewer]    |
| `.brix`            | File format       | Electron density map                | [Biostructure Viewer]    |
| `.mrc`             | File format       | Electron density volume                | [Biostructure Viewer]    |
| `.xplor`           | File format       | Electron density map                | [Biostructure Viewer]    |
| `.cub`, `.cube`    | File format       | 3D orbital or charge density grid | [Biostructure Viewer]    |
| `.dsn6 `           | File format       | Ribbons-format electron map                    | [Biostructure Viewer]    |

##### Simulation and topology formats

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.gro`             | File format       | GROMACS structure file                             | [Biostructure Viewer]    |
| `.top `            | File format       | GROMACS topology definition                              | [Biostructure Viewer]    |
| `.prmtop`          | File format       | AMBER force field topology                  | [Biostructure Viewer]    |
| `.parm7`           | File format       | AMBER topology file        | [Biostructure Viewer]    |
| `.psf `            | File format       | CHARMM or NAMD structure topology               | [Biostructure Viewer]    |
| `.cns`             | File format       | NMR structure description for CNS                | [Biostructure Viewer]    |

##### 3D mesh formats

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.obj `            | File format       | 3D mesh geometry                            | [Biostructure Viewer]    |
| `.ply`             | File format       | Polygon mesh data                            | [Biostructure Viewer]    |
| `.xyz`             | File format       | Atomic coordinates in XYZ format          | [Biostructure Viewer]    |


### Bioassays and experimental data

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------------------------|----------------------|
| `.xlsx`            | Data structure       | Assay plate data (auto-detected in file)                        | [Curves]             |

### Spectroscopy and NMR

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|--------------------------|
| `.jdx `            | File format       | Spectroscopy data (JCAMP-DX)     | [Chem Spectra Viewer]    |
| `.dx` | File format       | NMR spectroscopy data (JCAMP-DX variant)      | [NMRium]                 |
| `.nmrium`| File format       | NMR experiment data for NMRIum viewer     | [NMRium]                 |

### Biosignals

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------------------------|----------------------|
| `.edf`             | File format       | Multichannel biosignal recordings (e.g., EEG, EMG) | --                     |

### Clinical and regulatory data

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------------------------|----------------------|
| `.csv`             | Data structure       | SDTM-formatted datasets (auto-detected from folder structure)   | [Clinical Case]      |

### Geospatial formats

| Format <div style={{ width:90 }}></div> | Type <div style={{ width:150 }}></div> | Description <div style={{ width:400 }}></div> | Required plugins <div style={{ width:75 }}></div> |
|------------------|-------------------|------------------------------------------------|------------------|
| `.kml`       | File format       | Geospatial annotation file (Google Earth)               | [GIS]            |
| `.kmz`       | File format       | Compressed geospatial annotation file              | [GIS]            |
| `.geojson`         | File format       | Geospatial data encoded in JSON                      | [GIS]            |
| `.topojson`        | File format       | Compressed topological geospatial data         | [GIS]            |


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