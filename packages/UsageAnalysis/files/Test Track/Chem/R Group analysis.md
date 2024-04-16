1. Open smiles.csv (For example from Files -> Demo Files/chem/smiles.csv)
2. Go to **Chem > Analyze > R-groups analysis**. A dialog opens.
3. In the dialog, click the **MCS** button
4. Select the **Visual analysis** checkbox
5. Click the **OK** button.

**Expected**: balloon `Not enough R group columns to create trellis plot`

***

1. Open sar_small.csv
7. Go to **Chem > Analyze > R-groups analysis**. A dialog opens.
7. In the dialog, click the **MCS** button.
5. Click the **OK** button. **Expected**: a trellis plot with the results.
1. Run **R-groups analysis** once more.
1. Click **MCS**.
2. In the dialog **uncheck** the **Replace latest** checkbox. **Expected**: the second trellis plot should be displayed. In the grid, there should be two sets of columns resulting from the R-Group analysis.
1. Run **R-groups analysis** once more.
1. Click **MCS**.
2. In the dialog **check** the **Replace latest** checkbox. **Expected**: the latest results (columns and trellis plot) are replaced.
1. Run **R-groups analysis** once more.
1. At the bottom of the dialog, click the **Undo latest analysis**. **Expected**: the latest results (trellis and columns in the grid) should be removed. 
1. Run **R-groups analysis** without clicking the MCS. **Expected**: balloon
---
{
  "order": 2
}