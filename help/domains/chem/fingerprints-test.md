<!-- TITLE: Tests: Fingerprints -->
<!-- SUBTITLE: -->

# Tests: Fingerprints

Similarity measures, calculations that quantify the similarity of two molecules, and screening, 
a way of rapidly eliminating molecules as candidates in a substructure search, are both processes 
that use fingerprints. Fingerprints are a very abstract representation of certain structural features 
of a molecule.

## Testing scenario

1. Open *smiles.csv* file

1. Click on *"smiles"* column

1. Open *"Fingerprints"* dialog in *"Active"* tab on [Property Panel](../overview/navigation.md#properties) 
   (or use context menu for *"smiles"* column, **Chem | Fingerprintss...**)
   * *"Fingerprints"* dialog is open
   * Help switched to "Fingerprints" page

1. Execute *"Fingerprints"* dialog with default parameters
   * New column "Fingerprints" added to table
   * Column Fingerprints" contains numbers of bits and sets for each molecule

1. Open *"Fingerprints"* dialog again
   * History has saved parameters from previous run

1. Alternately execute dialog witch changing "Fingerprinter" field
   * When changing "Fingerprinter" field the dialog parameters are change
   * Required fields highlighted in red
   * If value in required field is no value, ```OK``` is not active for clicking
   * In result of dialog execution new column is created in table with values ​​corresponding to selected Fingerprinter

1. Test non-functional modules (UI, help, navigation, console)

See also:
 * [Molecular descriptor](../domains/chem/descriptors.md)
 * [Substructure Search Test](../tests/substructure-search-test.md)
 * [Get R Groups Test](../tests/get-r-groups-test.md)
