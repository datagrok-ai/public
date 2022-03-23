<!-- TITLE: Tests: R-group analysis test -->
<!-- SUBTITLE: -->

# Tests: R-group analysis

Gets R-groups of each molecule from input list in SMILES format around core.

## Testing scenarios

1. Open *smiles_mcs.csv* file
2. Open "R-groups Analysis" dialog from **Chem | R-groups Analysis**
3. Click on *"MCS"* button
    * In sketcher was drawn molecule structure, which corresponds to find MCS value

4. Execute "R-groups Analysis" dialog
    * Column "R0" is added to table, which corresponds to R-group "R0"

5. Count MCS separately (or use context menu for *"smiles"* column, **Chem | Find MCS**)

6. Copy calculated MCS value from console (**Tools | Console**) to sketcher text field in "R-groups Analysis" dialog
    * Molecule structure of described molecule smiles format was drawn in sketcher editor

7. Change value of "Column prefix" field to "test" and execute dialog
    * Column "test0" is added to table, which corresponds to R-group "R0"

See also:

* [Substructure search test](substructure-search-test.md)
* [Chem Sketcher test](chem-sketcher-test.md)
