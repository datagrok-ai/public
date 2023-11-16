# Monomer-Position viewer

Monomer-Position viewer is a matrix-like grid, where rows correspond to monomers and columns to positions. This viewer works in two modes: Mutation Cliffs and Invariant Map.

## Mutation Cliffs mode

Mutation Cliffs mode highlights the most signficant monomer-position pairs. The significance of the monomer-position is 
color-coded with blue for negative mean difference of the sequences containing monomer-position vs. rest of the sequences,
and red for positive mean difference. Color intensity is defined by p-value. The size of the circle represents the value of mean difference:
smaller circles represent low mean difference values; large circles represent hight mean difference values. Some of the 
cells contain number of unique sequences for given monomer-position, that form Mutation Cliff pairs.

## Invariant Map mode

Invariant Map mode shows number of sequences that contain monomer-position. Cells are color-coded based on the aggregation 
function on some numeric column. By default average scaled activity values for given monomer-position is used. Column and
aggregation function can be changed in viewer properties.

## Interactivity

Monomer-Position viewer selection follows default pattern as any other Datagrok viewer.

|                  |                                        |
|------------------|----------------------------------------|
| Click            | Single selection of Monomer-Position   |
| Ctrl+Click       | Toggle Monomer-Position selection state|
| Shift+Click      | Add Monomer-Position to selection      |
| Ctrl+Shift+Click | Deselect Monomer-Position              |

Hovering over any cell will show tooltip with distribution and statistics of the monomer-position and highlight corresponding
points on other viewers.
