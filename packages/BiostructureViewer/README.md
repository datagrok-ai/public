# BiostructureViewer

BiostructureViewer is a [package](https://datagrok.ai/help/develop/develop#packages) for
the [Datagrok](https://datagrok.ai) platform for biological structures visualization.

## Formats

The following file formats are supported:

* Molecular structures (mmCIF, PDB, PQR, GRO, SDF, MOL2, MMTF)
* Density volumes (MRC/MAP/CCP4, DX/DXBIN, CUBE, BRIX/DSN6, XPLOR/CNS)

## Semantic Types

The BiostructureViewer package enables detecting and handling the _Molecule3D_ semantic type.

## Cell renderers

Grid cells of the columns of _Molecule3D_ semantic type are subjected to draw with a designated NGL based cell renderer.
Mouse click on a cell is handled opening a [BiostructureViewer](#BiostructureViewer) to explore the structure details.

![pdb_data](../../help/uploads/gifs/biostructure-viewer-pdb-data.gif)

## Viewers

Both [BiostructureViewer](#BiostructureViewer) and [NglViewer](#NglViewer) utilize the same way of obtaining identical a
structural data into the viewer from different sources (in order of priority) controlled by
properties of the 'Data' category:

1. A PDB string value of the 'pdb' property.
2. A PDB string value of the data frame tag named of 'pdbTag' property.
3. A PDB string value of the data frame '.pdb' tag in case 'pdbTag' property is not specified.

### <a  name="NglViewer">NglViewer</a>

[NglViewer](./src/viewers/ngl-viewer.ts) is a Datagrok [DG.JsViewer](../../js-api/src/viewer.ts) derived
component based on [NGL Viewer](https://nglviewer.org/) library developed by Alexander Rose.

![ngl-viewer](../../help/uploads/gifs/ngl-viewer-open-PDB.gif)

### <a name="BiostructureViewer">BiostructureViewer</a>

[BiostructureViewer](./src/viewers/molstar-viewer.ts) is a Datagrok [DG.JsViewer](../../js-api/src/viewer.ts) derived
component based on [RCSB PDB implementation](https://github.com/molstar/rcsb-molstar) of
[Mol*](https://github.com/molstar/molstar).
Documentation of the Mol* project can be found [here](https://molstar.org/docs/).

Exposed 'Style' category properties allow customizing the viewer appearance representation as
cartoon, backbone, ball+stick, licorice, hyperball, surface.

The structural data can be put into the viewer from different sources (in order of priority) controlled by properties
of the 'Data' category:

1. A PDB string value of the 'pdb' property.
2. A PDB string value of the data frame tag named of 'pdbTag' property.
3. A PDB string value of the data frame '.pdb' tag in case 'pdbTag' property is not specified.

#### References

1. David Sehnal, Sebastian Bittrich, Mandar Deshpande, Radka Svobodová, Karel Berka,
   Václav Bazgier, Sameer Velankar, Stephen K Burley, Jaroslav Koča, Alexander S Rose:
   Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures,
   Nucleic Acids Research, 2021; [doi.org/10.1093/nar/gkab314](https://doi.org/10.1093/nar/gkab314).

### NglViewer

Exposed 'Style' category properties allow customizing the viewer appearance representation as
cartoon, backbone, ball+stick, licorice, hyperball, surface.

## File viewers

Files with handled extension are previewed and open with the [BiostructureViewer](#BiostructureViewer).

This package implements [custom file viewers](../../help/develop/how-to/create-custom-file-viewers.md)
for all supported file types. This is how it looks in action:

![viewers](../../help/access/file-shares-file-viewers.gif)

See also:

* [NGL Viewer](https://nglviewer.org/)
* [Scripting](https://datagrok.ai/help/compute/scripting)
* [Files](https://datagrok.ai/help/access/connectors/files)
* [File Shares](https://datagrok.ai/help/access/file-shares)
* [RCSB PDB implementation](https://github.com/molstar/rcsb-molstar)

## Versions

@rcsb/rcsb-molstar v1.8.7 is the last dependent on rxjs of version 6.x.x
and it has dependency on "[molstar](https://github.com/molstar/molstar)": "^2.4.1"
