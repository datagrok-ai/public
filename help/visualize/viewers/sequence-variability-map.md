# Sequence Variability Map viewer

The Sequence Variability Map Viewer presents data in a matrix-like grid where rows represent monomers and columns represent positions. This viewer operates in two distinct modes: Mutation Cliffs and Invariant Map.

## Mutation Cliffs mode

Mutation Cliffs Mode highlights significant monomer-position pairs that differ at only one specific position, with one sequence containing the row-corresponding monomer at that position. The significance of each monomer-position pair is visually encoded through color:

* **Blue shades**: Indicate negative mean activity difference between sequences containing the monomer-position versus all other sequences

* **Red shades**: Indicate positive mean activity difference between sequences containing the monomer-position versus all other sequences

The size of each circle represents the number of unique sequences associated with a particular mutation cliff.

When you select a mutation cliff, the viewer automatically selects all sequence pairs that differ solely at that position and monomer.

You can further inspect the mutation cliffs from the context panel, which displays mutations and their corresponding sequences.

## Invariant Map mode

Invariant Map Mode displays the number of sequences containing each specific monomer-position combination. Cells are color-coded based on an aggregation function applied to a numeric column (with activity average as the default metric).

Selecting a monomer-position combination in this mode will select all sequences containing the specified monomer at the given position.

## Customization

You can modify the Sequence Variability Map Viewer by changing the following properties from the property panel:

- **Sequence**: The column containing the sequence data.
- **Activity**: The column containing the activity data.
- **Activity Scaling**: The scaling method for the activity data.
- **Target**: The target column to be used as internal filter.
- **Aggregations**: Additional columns to aggregate in the table.

    ### Mutation Cliffs
    
    - **Min Activity Delta**: The minimum activity difference to be considered in mutation cliffs (0 by default).
    - **Max Mutations**: The maximum number of point mutations to be considered in mutation cliffs (1 by default).

    ### Invariant Map

    - **Value**: Activity column to be used for invarian map.
    - **Value Aggregation**: Aggregation function to be applied to the activity column.
    - **Color**: Color column to be used for invariant map, defaults to activity column.
    - **Color Aggregation**: Aggregation function to be applied to the color column.
    - **Lower/Mid/Upper Bound Color**: Bounds for the color scale.
    - **Log Scale Color**: Use logarithmic scale for color.
    - **Custom Color Range**: Custom color range for the color scale (instead of column min/max).

## Interactivity

The Sequence Variability Map Viewer follows Datagrok's standard selection patterns:

|                  |                                        |
|------------------|----------------------------------------|
| Click            | Single selection of Monomer-Position   |
| Ctrl+Click       | Toggle Monomer-Position selection state|
| Shift+Click      | Add Monomer-Position to selection      |
| Ctrl+Shift+Click | Deselect Monomer-Position              |

Hovering over any populated cell displays a tooltip showing distribution statistics for the monomer-position combination and highlights the corresponding rows in other active viewers.

[Sequence Variability Map](./img/SVM-viewer.gif)