# Supported formats

## Tabular formats

| Extension  | Description                                                    | Source                 |
|------------|----------------------------------------------------------------|------------------------|
| .txt       | Plain text                                                     | Core                   |
| .csv       | Comma-separated values                                         | Core                   |
| .tsv       | Tab-separated values                                           | Core                   |
| .xml       | Extensible Markup Language                                     | Core                   |
| .json      | JavaScript Object Notation                                     | Core                   |
| .HTML      | HyperText Markup Language                                      | Core                   |
| .xlsx      | Excel                                                          | Core                   |
| .edf       | European Data Format                                           | Core                   |
| .sas7bdat  | SAS database file                                              | Core                   |
| .kxl       | KeyCreator eXtensible Language                                 | Core   <!--check!!!--> |
| .rds, .rda | R Data Format                                                  | Core                   |
| .h5        | Hierarchical Data Format                                       | Core                   |
| .nc        | NetCDF                                                         | Core                   |
| .mat       | MATLAB MAT                                                     | Core                   |
| .d42       | Datagrok [project](../../datagrok/concepts/project/project.md) | Core                   |
| .zip       | ZIP archive (for supported types)                              | Core                   |
| .gz, .gzip | gzip                                                           | Core                   |
| .tar       | Tape archive                                                   | Core                   |
| netCDF     | network Common Data Form                                       | Core  <!--check!!!-->  |
| .kml, .kmz | Keyhole Markup Language (geographic annotations)               | [GIS]                  |
| .parquet   | Efficient data storage and retrieval                           | [Arrow]                |
| .feather   | portable file format for storing Arrow tables                  | [Arrow]                |
| .sqlite    | SQLite database                                                | [Dendrogram]           |

## Molecule structure formats

| Extension | Description                                                   | Source |
|-----------|---------------------------------------------------------------|--------|
| .smi    | Structure-data file                                           | [Chem] |
| .sdf      | Structure-data file                                           | [Chem] |
| .mol      | MDL Molfile                                                   | [Chem] |
| .mol2     | SYBYL molecule representation                                 | [Chem] |

## Molecule3D structure formats

| Extension | Description                                                   | Source               | Viewer         |
|-----------|---------------------------------------------------------------|----------------------|----------------|
| .brix     | (preview) Electron Density Map                                | [BiostructureViewer] |                |
| .ccp4     | Collaborative Computational Project Number 4 electron density | [BiostructureViewer] | [NGL]          |
| .cif      | (preview) Crystallographic Information File                   | [BiostructureViewer] |                |
| .cifCore  | Crystallographic Information Core File                        | [BiostructureViewer] | [Biostructure] |
| .cns      | Density ?                                                     | [BiostructureViewer] | [NGL]          |
| .cub      | (preview) Orbital/Density Values on a 3D Grid                 | [BiostructureViewer] |                |
| .cube     | (preview) Orbital/Density Values on a 3D Grid                 | [BiostructureViewer] |                |
| .dsn6     | (preview) Ribbons FRODO Map                                   | [BiostructureViewer] |                |
| .dx       | (preview) Chemical Spectroscopy Format                        | [BiostructureViewer] |                |
| .dxbin    | (preview) .dxbin                                              | [BiostructureViewer] |                |
| .ent      | (preview) Brookhaven PDB Molecule                             | [BiostructureViewer] |                |
| .gro      | GROMACS                                                       | [BiostructureViewer] | [Biostructure] |
| .mcif     | (preview) Macromolecular Crystallographic Information File    | [BiostructureViewer] |                |
| .mmCif    | Macromolecular Crystallographic Information File              | [BiostructureViewer] | [Biostructure] |
| .mmtf     | Macromolecular Transmission Format                            | [BiostructureViewer] | [NGL]          |
| .mrc      | (preview) Electron Density Map                                | [BiostructureViewer] |                |
| .obj      | Wavefront .obj file                                           | [BiostructureViewer] | [NGL]          |
| .parm7    | (preview) AMBER Parameter/Topology File                       | [BiostructureViewer] |                |
| .pdb      | Protein Data Bank                                             | [BiostructureViewer] | [Biostructure] |
| .pdbqt    | Protein Data Bank, Partial Charge (Q), & Atom Type (T)        | [BiostructureViewer] | [Biostructure] |
| .ply      | Polygon file format                                           | [BiostructureViewer] | [NGL]          |
| .psf      | (preview) Protein Structure File                              | [BiostructureViewer] |                |
| .sd       | (preview) Structure-Data                                      | [BiostructureViewer] |                |
| .prmtop   | Parameter / Topology (AMBER)                                  | [BiostructureViewer] | [NGL]          |
| .top      | Gromacs Topology File                                         | [BiostructureViewer] | [NGL]          |
| .xplor    | (preview) Electron Density Map                                | [BiostructureViewer] |                |
| .xyz      | (preview) Cartesian Coordinates Format                        | [BiostructureViewer] |                |

## Macromolecule formats

| Extension | Description                                    | Source |
|-----------|------------------------------------------------|--------|
| .fasta    | FASTA format (fa, fna, ffn, faa, frn, fa, fst) | [Bio]  |

Note: The platform detector also recognizes semantic types of Macromolecule in files that are in [Tabular format](#tabular-formats).

## Miscellaneous

| Extension      | Description                        | Source              |
|----------------|------------------------------------|---------------------|
| .ipynb         | Jupyter Notebook                   | [Notebooks]         |
| .ipv           | Interactive differential equations | [Dendrogram]        |
| .pdf           | Portable document format           | [FileEditors]       |
| .docx          | MS Word                            | [FileEditors]       |
| .jdx,dx,nmrium | NMR                                | [NMRium]            |
| .jdx           | Spectroscopy                       | [ChemSpectraViewer] |


See also:

[Chem]: <https://github.com/datagrok-ai/public/tree/master/packages/Chem#readme>
[BiostructureViewer]: <https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer#readme>
[NGL]: ../../visualize/viewers/ngl.md
[Biostructure]: ../../visualize/viewers/biostructure.md

[GIS]: https://github.com/datagrok-ai/public/tree/master/packages/GIS#readme
[Arrow]: https://github.com/datagrok-ai/public/tree/master/packages/Arrow#readme
[Dendrogram]: https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram#readme
[FileEditors]: https://github.com/datagrok-ai/public/tree/master/packages/FileEditors#readme
[DiffStudio]: https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio#readme
[SQLite]: https://github.com/datagrok-ai/public/tree/master/packages/SQLite#readme
[Notebooks]: https://github.com/datagrok-ai/public/tree/master/packages/Notebooks#readme
[NMRium]: https://github.com/datagrok-ai/public/tree/master/packages/NMRium#readme
[ChemSpectraViewer]: https://github.com/datagrok-ai/chem-spectra-viewer

[//]: # ([Notebooks]: https://github.com/datagrok-ai/public/tree/master/packages/Notebooks#readme)
