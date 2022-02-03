<!-- TITLE: Tests: Similarity analysis -->
<!-- SUBTITLE: -->

# Testing scenarios

1. Open *smiles.csv* file

2. Open "Similarity analysis" dialog from **Chem | Similarity Analysis**

   * Dialog shows number of not not represent valid molecules
   * You can choose actions with not represent valid molecules (Filter, Select or Delete)

3. Click on *"Select* button

   * Rows with not represent valid molecules are selected

4. Click on *"Filter* button

   * Rows with not represent valid molecules are filtered

5. Click on *"Delete* button

   * Rows with not represent valid molecules are deleted from table

6. Click on *"Re-Run* button

   * Similarity analysis is performed
   * Columns "x" and "y" with coordinates are added to table
   * [Scatter Plot](../../visualize/viewers/scatter-plot.md) added to table layout
   * Special column "~ RDKFingerprint" added to table

See also:

* [Substructure search test](substructure-search-test.md)
* [Chem Sketcher Test](chem-sketcher-test.md)
