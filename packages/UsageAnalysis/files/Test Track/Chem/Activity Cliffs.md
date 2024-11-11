1. Open:
   * Linked datasets (should be tested on smiles (smiles), molV2000 (mol1K) and molV3000 (ApprovedDrugs2015) formats)
   * smiles_2_columns.csv (TODO: add to linked datasets)
2. Run **Chem > Analyze > Activity cliffs**. A dialog opens.
3. Click OK to run a function with default parameters
4. Check 'Show only cliffs' - only points with cliffs should be left on scatter plot
5. Click on the link with number of cliffs - grid with cliffs should appear on the bottom
6. Click on the first row of the grid - the scatter plot is zoomed to the cliff, property panel is shown with the pair of molecules (highlighting uncommon parts), activity difference. When hovering over line the tooltip with the same information should be showed.
7. Double click to unzoom the scatter plot. And zoom any cluster to see any line. Click on the line - the scatter plot will be zoomed. And corresponding line will become current in the grid below.
9. Run **Activity cliffs** one more time.
10. Change the parameters arbitrarily
11. Click OK to run a function with edited parameters.
---
{
  "order": 7,
  "datasets": [
    "System:DemoFiles/SPGI.csv",
    "System:DemoFiles/chem/smiles.csv",
    "System:AppData/Chem/mol1K.sdf",
    "System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf",
    "System:AppData/UsageAnalysis/test_datasets/smiles_2_columns.csv"
  ]
}
