<!-- TITLE: Tests: R-Group Analysis Test -->
<!-- SUBTITLE: -->

# Tests: R-Group Analysis

Gets R-groups of each molecule from input list in SMILES format around core.

## Testing scenarios

1. Open *smiles_mcs.csv* file 

1. Open "R-groups Analysis" dialog from **Chem | R-groups Analysis**

1. Click on *"MCS"* button 
   * In sketcher was drawn molecule structure, which corresponds to find MCS value
   
1. Execute "R-groups Analysis" dialog
   * Column "R0" is added to table, which corresponds to R-group "R0"
   
1. Count MCS separately (or use context menu for *"smiles"* column, **Chem | Find MCS**)  
   
1. Copy calculated MCS value from console (**Tools | Console**) to sketcher text field in "R-groups Analysis" dialog
   * Molecule structure of described molecule smiles format was drawn in sketcher editor
   
1. Change value of "Column prefix" field to "test" and execute dialog
      * Column "test0" is added to table, which corresponds to R-group "R0"

See also:
 * [Substructure Search Test](../tests/substructure-search-test.md)
 * [Chem Descriptors Test](../tests/chem-descriptors-test.md)
 * [Chem Sketcher Test](../tests/chem-scetcher-test.md)