# Supported formats

Datagrok supports a wide range of file formats across domains such as cheminformatics, bioinformatics, clinical, and geospatial data. It also handles a set of molecular and macromolecular notations.

In addition to detecting formats by their extension, Datagrok uses content-driven data handling to choose the most effective way to display the data.

## Structured and analytical data

### Tabular and semi-structured data

| Extension         | Description                                    | Source    |
|------------------|------------------------------------------------|--------------------|
| .csv             | Comma-separated values                         | Core               |
| .tsv             | Tab-separated values                           | Core               |
| .txt             | Plain text                                     | Core               |
| .xlsx            | Microsoft Excel                                | Core               |
| .d42             | Datagrok [project]                             | Core               |
| .json            | JavaScript Object Notation                     | Core               |
| .xml             | Extensible Markup Language                     | Core               |
| .html            | HyperText Markup Language                      | Core               |

### Statistical and scientific formats

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .rds, .rda        | R serialized data                              | Core               |
| .mat              | MATLAB data format                             | Core               |
| .h5               | Hierarchical Data Format                       | Core               |
| .nc, .netcdf      | Network Common Data Form                       | Core               |
| .sas7bdat         | SAS database file                              | Core               |
| .kxl              | KeyCreator eXtensible Language                 | Core               |

### Columnar and database formats

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .parquet          | Columnar storage format                        | [Arrow]              |
| .feather          | Arrow-based binary format                      | [Arrow]              |
| .sqlite           | SQLite database                                | [SQLite]         |

### Archives and packaging

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .zip              |ZIP archive (auto-extracted for supported types)| Core               |
| .tar              | Tape archive                                   | Core               |
| .gz, .gzip        | Gzip-compressed files                          | Core               |


## Notebooks and documents

### Notebooks and computations

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .ipynb            | Jupyter Notebook                               | [Notebooks]        |
| .ivp              | Interactive differential equations             | [Diff Studio]      |

### Docs and reports

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .pdf              | Portable Document Format                       | [File Editors]     |
| .docx             | Microsoft Word Document                        | [File Editors]     |


## Domain-specific formats

### Molecular structures

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .sdf, .sd         | Structure-data file                            | [Chem]             |
| .mol              | MDL Molfile                                    | [Chem]             |
| .mol2             | SYBYL molecule representation                  | [Chem]             |
| .smi              | SMILES in file format                          | [Chem]             |

### Molecular and macromolecular notations (used in strings or columns)

| Notation          | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| SMILES            | Simplified Molecular Input Line Entry System   | [Chem]             |
| SMARTS            | SMiles ARbitrary Target Specification          | [Chem]             |
| InChI, InChIKey   | IUPAC chemical identifiers                     | [Chem]             |
| HELM              |Hierarchical Editing Language for Macromolecules| [Bio]       |
| SMIRKS            | Reaction transformations                       | [Chem]             |
| BILN              | Biopolymer Line Notation                       | [Bio]      |

### Macromolecules and sequences

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .fasta            | FASTA format (fa, fna, ffn, faa, frn, fa, fst) | [Bio]              |

### Phylogenetic trees

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|-----------|
| .nwk, .newick     | Newick format for phylogenetic trees           | [PhyloTree Viewer]     |

### 3D structures and density maps

#### Macromolecular structure formats

| Extension         | Description                                    | Source    |
|------------------|------------------------------------------------|--------------------|
| .pdb              | Protein Data Bank                              | [Biostructure Viewer] |
| .pdbqt     | Protein Data Bank, Partial Charge (Q), & Atom Type (T)| [Biostructure Viewer] |
| .mcif |(preview) Macromolecular Crystallographic Information  File | [Biostructure Viewer] |
| .mmcif          | Macromolecular Crystallographic Information  File| [Biostructure Viewer] |
| .ent              | (preview) Brookhaven PDB Molecule              | [Biostructure Viewer] |
| .cif              | Crystallographic Information File              | [Biostructure Viewer] | 
| .cifCore          | Crystallographic Information Core File         | [Biostructure Viewer] |
| .mmtf             | Macromolecular Transmission Format             | [Biostructure Viewer] |

#### Electron density and volumetric data

| Extension         | Description                                    | Source    |
|------------------|------------------------------------------------|--------------------|
| .ccp4            | Collaborative Computational Project Number 4  | [Biostructure Viewer] |
| .brix            | (preview) Electron Density Map                | [Biostructure Viewer] |
| .mrc             | (preview) Electron Density Map                | [Biostructure Viewer] |
| .cub, .cube      | (preview) Orbital/Density Values on a 3D Grid | [Biostructure Viewer] |
| .dx, .dxbin      | (preview) Chemical spectroscopy format        | [Biostructure Viewer] |
| .dsn6            | (preview) Ribbons FRODO map                   | [Biostructure Viewer] |
| .xplor           | (preview) Electron density map                | [Biostructure Viewer] |

#### Simulation and topology formats

| Extension         | Description                                    | Source    |
|------------------|------------------------------------------------|--------------------|
| .gro          | GROMACS structure                            | [Biostructure Viewer]   |
| .top          | GROMACS topology                             | [Biostructure Viewer]   |
| .prmtop       | AMBER topology                               | [Biostructure Viewer]   |
| .parm7        | (preview) AMBER Parameter/Topology File      | [Biostructure Viewer]   |
| .psf          | (preview) Protein Structure File             | [Biostructure Viewer]   |
| .cns          | Crystallography and NMR System               | [Biostructure Viewer]   |

#### 3D Mesh Formats

| Extension     | Description                                  | Source                |
|---------------|----------------------------------------------|------------------------|
| .obj          | Wavefront 3D object                          | [Biostructure Viewer]     |
| .ply          | Polygon file format                          | [Biostructure Viewer]     |
| .xyz          | (preview) Cartesian Coordinate Format        | [Biostructure Viewer]     |

### Spectroscopy and NMR

| Extension         | Description                                    | Source    |
|------------------|------------------------------------------------|----------------------|
| .jdx              | Spectroscopy                                  | [Chem Spectra Viewer]|
| .jdx, .dx, .nmrium| NMR                                           | [NMRium]             |

### Clinical and regulatory data

* **SDTM** (Study Data Tabulation Model) datasets: the [Clinical Case] application autodetects SDTM-formatted datasets in folders containing `.csv` files and loads them for analysis.

### Assay plate data

* [Curves] autodetects layout, concentration, and activity layers in `.xlsx` files, and visualizes plate content as fitted dose-response curves.

### Biomedical signals

| Extension         | Description                                    | Source    |
|-------------------|------------------------------------------------|--------------------|
| .edf              | European Data Format for biosignals            | Core               |

### Geospatial formats

| Extension         | Description                                    | Source    |
|------------------|------------------------------------------------|--------------------|
| .kml, .kmz        | Keyhole Markup Language                       | [GIS]                |
| .geojson          | GeoJSON (geospatial data)                     | [GIS]                |
| .topojson         | TopoJSON (topological geospatial data)        | [GIS]                |

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