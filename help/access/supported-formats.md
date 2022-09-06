<!-- TITLE: File -->
<!-- SUBTITLE: -->

# Supported formats

Datagrok supports the following formats.

## Tabular formats

| Extension  | Description                                      | Source |
|------------|--------------------------------------------------|--------|
| .txt       | Plain text                                       | Core   |
| .csv       | Comma-separated values                           | Core   |
| .tsv       | Tab-separated values                             | Core   |
| .xml       | Extensible Markup Language                       | Core   |
| .json      | JavaScript Object Notation                       | Core   |
| .HTML      | HyperText Markup Language                        | Core   |
| .xlsx      | Excel                                            | Core   |
| .edf       | European Data Format                             | Core   |
| .sas7bdat  | SAS database file                                | Core   |
| .kml, .kmz | Keyhole Markup Language (geographic annotations) | Core   |
| .kxl       | KeyCreator eXtensible Language                   | Core   |<!--check!!!-->
| .rds, .rda | R Data Format                                    | Core   |
| .h5        | Hierarchical Data Format                         | Core   |
| .nc        | NetCDF                                           | Core   |
| .mat       | MATLAB MAT                                       | Core   |
| .d42       | Datagrok [project](../datagrok/project.md)       | Core   |
| .zip       | ZIP archive (for supported types)                | Core   |
| .gz, .gzip | gzip                                             | Core   |
| .tar       | Tape archive                                     | Core   |
| .ipynb     | Jupyter Notebook                                 | Core   |
| netCDF     | network Common Data Form                         | Core   | <!--check!!!-->

![File browsing](files-browser.gif "Files Browser") <!--rename-->

## Molecular structure formats

| Extension | Description                        | Source             |
|-----------|------------------------------------|--------------------|
| .cif      | Crystallographic Information File  | [NglViewer] plugin |
| .pdb      | Protein Data Bank                  | [NglViewer] plugin |
| .pqr      | PQR                                | [NglViewer] plugin |
| .gro      | GROMACS                            | [NglViewer] plugin |
| .sdf      | Structure-data file                | [NglViewer] plugin |
| .mol      | MDL Molfile                        | [NglViewer] plugin |
| .mol2     | SYBYL molecule representation      | [NglViewer] plugin |
| .mmtf     | Macromolecular Transmission Format | [NglViewer] plugin |

See also:

* [File shares](connect-a-file-share.md)

[NglViewer]: https://github.com/datagrok-ai/public/tree/master/packages/NglViewer#readme

[//]: # ([Notebooks]: https://github.com/datagrok-ai/public/tree/master/packages/Notebooks#readme)
