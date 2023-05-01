# <a name="NglViewerl">NglViewer</a>

[NglViewer](../../../packages/BiostructureViewer/src/viewers/ngl-viewer.ts) is
a Datagrok [DG.JsViewer](../../js-api/src/viewer.ts) derived component based on
the NGL Viewer [nglviewer.org](https://nglviewer.org/) library developed by Alexander Rose.

NglViewer obtains a structural data into the viewer from different sources
(in order of priority) controlled by properties of the 'Data' category:

1. A PDB string value of the 'pdb' property.
2. A PDB string value of the data frame tag named of 'pdbTag' property.
3. A PDB string value of the data frame '.pdb' tag in case 'pdbTag' property is not specified.

The 'Behaviour' category properties allow control interactivity with dataframe
selected, current, and mouse over rows, showing the ligands along with a biostructure.

## <a name="">Ligands</a>

The Datagrok tool is used to explore docking results. The NglViewer feature can load a structure from a file
using the 'open...' button, or it can use structure data from the underlying data frame tag, '.pdb' by default.

The NglViewer detects a column of semantic type 'Molecule' or specified with a property 'ligandColumnName' and
displays its data along with the structure handling data frame rows according to the 'Behaviour' category properties.

![biostructure-viewer-ngl-viewer-ligand-interactivity](../../uploads/gifs/ngl-viewer-ligand-interactivity.gif)

### Behaviour Category Properties to show Ligands

| Property Name           | Default | Description                                                                           |
|-------------------------|---------|---------------------------------------------------------------------------------------|
| ligandColumnName        | string  | The name of the ligand source column <br/> (autodetected by semantic type 'Molecule') |
| showSelectedRowsLigands | false   | Controls the display of ligands of selected rows                                      |
| showCurrentRowLigand    | true    | Controls the display of the ligand of the current row                                 |
| showMouseOverRowLigand  | true    | Controls the display of the ligand of a mouse-over row                                |

If only one ligand source is checked, then the ligand is displayed with a full-color ball+stick
representation, and selected rows ligands with colors from a color scheme scale. If there are multiple ligand sources,
then the ligand of the current row is displayed in 'green', the mouse-over row ligand is displayed in 'light gray',
and selected rows ligands are displayed in orange.
