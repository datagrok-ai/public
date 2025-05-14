# Supported formats

Datagrok supports over 50 file formats, including common and domain-specific types across cheminformatics, bioinformatics, clinical data, and more, with a primary focus on structured tabular data analysis.
In addition, Datagrok supports domain-specific notations and understands data formatted in specific ways, e.g. assay plates and SDTM datasets, enabling content-aware parsing and visualization.

## Structured and analytical data

### Tabular and semi-structured data

| Extension         | Description                                   | 
|------------------|------------------------------------------------|
| .csv             | Comma-separated values                         | 
| .tsv             | Tab-separated values                           |
| .txt             | Plain text                                     |
| .xlsx            | Microsoft Excel                                |
| .d42             | Datagrok [project]                             |
| .json            | JavaScript Object Notation                     |
| .xml             | Extensible Markup Language                     |
| .HTML            | HyperText Markup Language                      |

### Statistical and scientific formats

| Extension         | Description                                    | 
|-------------------|------------------------------------------------|
| .rds, .rda        | R serialized data                              |
| .mat              | MATLAB data format                             | 
| .h5               | Hierarchical Data Format                       |  
| .nc, .netcdf      | Network Common Data Form                       | 
| .sas7bdat         | SAS database file                              |
| .kxl              | KeyCreator eXtensible Language                 | 

### Columnar and database formats

| Extension         | Description                                    | Requirements    |
|-------------------|------------------------------------------------|-----------------|
| .parquet          | Columnar storage format                        | [Arrow]         |
| .feather          | Arrow-based binary format                      | [Arrow]         |
| .sqlite           | SQLite database                                | [SQLite]        |

### Archives and packaging

| Extension         | Description                                    |
|-------------------|------------------------------------------------|
| .zip              |ZIP archive (auto-extracted for supported types)| 
| .tar              | Tape archive                                   |  
| .gz, .gzip        | Gzip-compressed files                          |


## Notebooks and documents

### Notebooks and computations

| Extension         | Description                                    | Requirements    |
|-------------------|------------------------------------------------|-----------------|
| .ipynb            | Jupyter Notebook                               | [Notebooks]     |
| .ivp              | Interactive differential equations             | [Diff Studio]   |

### Docs and reports

| Extension         | Description                                    | Requirements    |
|-------------------|------------------------------------------------|-----------------|
| .pdf              | Portable Document Format                       | [File Editors]  |
| .docx             | Microsoft Word Document                        | [File Editors]  |


## Domain-specific formats

### Cheminformatics

#### Molecular structures

| Extension         | Description                                    | Requirements       |
|-------------------|------------------------------------------------|--------------------|
| .sdf, .sd         | Structure-data file                            | [Chem]             |
| .mol              | MDL Molfile                                    | [Chem]             |
| .mol2             | SYBYL molecule representation                  | [Chem]             |
| .smi              | Structure-data file                            | [Chem]             |

### Bioinformatics

#### Macromolecules and sequences

| Extension         | Description                                    | Requirements       |
|-------------------|------------------------------------------------|--------------------|
| .fasta            | FASTA format (fa, fna, ffn, faa, frn, fa, fst) | [Bio]              |

#### Phylogenetic trees

| Extension         | Description                                    | Requirements      |
|-------------------|------------------------------------------------|-------------------|
| .nwk, .newick     | Newick format for phylogenetic trees           | [PhyloTree Viewer]|

#### 3D structures and density maps

##### Macromolecular structure formats

| Extension         | Description                                    | Requirements    |
|------------------|------------------------------------------------|--------------------|
| .pdb              | Protein Data Bank                              | [Biostructure Viewer] |
| .pdbqt     | Protein Data Bank, Partial Charge (Q), & Atom Type (T)| [Biostructure Viewer] |
| .cif    | (preview) Crystallographic Information File              | [Biostructure Viewer] | 
| .cifCore          | Crystallographic Information Core File         | [Biostructure Viewer] |
| .mcif |(preview) Macromolecular Crystallographic Information  File | [Biostructure Viewer] |
| .mmcif          | Macromolecular Crystallographic Information  File| [Biostructure Viewer] |
| .mmtf             | Macromolecular Transmission Format             | [Biostructure Viewer] |
| .ent              | (preview) Brookhaven PDB Molecule              | [Biostructure Viewer] |

##### Electron density and volumetric data

| Extension         | Description                                    | Requirements    |
|------------------|------------------------------------------------|--------------------|
| .ccp4            | Collaborative Computational Project Number 4  | [Biostructure Viewer] |
| .brix            | (preview) Electron Density Map                | [Biostructure Viewer] |
| .mrc             | (preview) Electron Density Map                | [Biostructure Viewer] |
| .xplor           | (preview) Electron density map                | [Biostructure Viewer] |
| .cub, .cube      | (preview) Orbital/Density Values on a 3D Grid | [Biostructure Viewer] |
| .dsn6            | (preview) Ribbons FRODO map                   | [Biostructure Viewer] |
<!--| .dx, .dxbin      | (preview) Chemical spectroscopy format        | [Biostructure Viewer] |-->

##### Simulation and topology formats

| Extension         | Description                                    | Requirements    |
|------------------|------------------------------------------------|--------------------|
| .gro          | GROMACS structure                            | [Biostructure Viewer]   |
| .top          | GROMACS topology                             | [Biostructure Viewer]   |
| .prmtop       | AMBER Parameter/Topology File                | [Biostructure Viewer]   |
| .parm7        | (preview) AMBER Parameter/Topology File      | [Biostructure Viewer]   |
| .psf          | (preview) Protein Structure File             | [Biostructure Viewer]   |
| .cns          | Crystallography and NMR System               | [Biostructure Viewer]   |

##### 3D mesh formats

| Extension     | Description                                  | Requirements                |
|---------------|----------------------------------------------|------------------------|
| .obj          | Wavefront 3D object                          | [Biostructure Viewer]     |
| .ply          | Polygon file format                          | [Biostructure Viewer]     |
| .xyz          | (preview) Cartesian Coordinate Format        | [Biostructure Viewer]     |

### Spectroscopy and NMR

| Extension         | Description                                    | Requirements    |
|------------------|------------------------------------------------|----------------------|
| .jdx              | Spectroscopy                                  | [Chem Spectra Viewer]|
| .jdx, .dx, .nmrium| NMR                                           | [NMRium]             |

### Clinical and regulatory data

| Extension   | Description                                                          | Requirements |
|-------------|----------------------------------------------------------------------|--------------|
| .edf        | European Data Format for biosignals                                  |              |
| .csv        | SDTM-formatted datasets (auto-detected from folder structure)        | [Clinical Case] |
| .xlsx       | Assay plate data (auto-detected in file)                             | [Curves]       |

### Geospatial formats

| Extension         | Description                                    | Requirements    |
|------------------|------------------------------------------------|--------------------|
| .kml, .kmz        | Keyhole Markup Language                       | [GIS]                |
| .geojson          | GeoJSON (geospatial data)                     | [GIS]                |
| .topojson         | TopoJSON (topological geospatial data)        | [GIS]                |

## Molecular and macromolecular notations (used in strings or columns)

| Notation          | Description                                    | Requirements    |
|-------------------|------------------------------------------------|--------------------|
| SMILES            | Simplified Molecular Input Line Entry System   | [Chem]             |
| SMARTS            | SMiles ARbitrary Target Specification          | [Chem]             |
| InChI, InChIKey   | IUPAC chemical identifiers                     | [Chem]             |
| HELM              |Hierarchical Editing Language for Macromolecules| [Bio]       |
| SMIRKS            | Reaction transformations                       | [Chem]             |
| BILN              | Biopolymer Line Notation                       | [Bio]      |

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