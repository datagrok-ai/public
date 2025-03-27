# Most Potent Residues

The Most Potent Residues viewer identifies the most potent monomer at each position and displays relevant statistics. Statistics include mean difference, p-value, count and ratio. The algorithm analyzes activity distributions with position invariants to identify monomers with the highest activity for each position. You can adjust the viewer through the property panel to target either high or low activity values.

Numbers in each cell show the quantity of sequences with the given monomer at that position. Circle size represents the mean difference between sequences with the given monomer at that position and all other sequences. Circle color indicates the p-value of the mean difference.

## Interactivity

Selecting a monomer-position will select all sequences containing that specific monomer at the given position.

The Most Potent Residues viewer follows Datagrok's standard selection pattern:

|                  |                                        |
|------------------|----------------------------------------|
| Click            | Single selection of Monomer-Position   |
| Ctrl+Click       | Toggle Monomer-Position selection state|
| Shift+Click      | Add Monomer-Position to selection      |
| Ctrl+Shift+Click | Deselect Monomer-Position              |

Hovering over a monomer-position cell displays a tooltip with distribution and statistics of the monomer-position and highlights corresponding points on other viewers.

[Most Potent Residues](./img/MPR.gif)