1. The sample_HELM.csv file is too large, so extract a subset of 50 rows
2. Add a new column with a formula RandBetween(0,5) - this will be the clusters column.
3. On the menu ribbon, open **Bio** > **MSA**.
4. Set Cluster to the new column with a formula (from the step 2).
5. Check that **Alignment parameters** button adds input parameters to the dialog properly   
6. Check the new column.

* Everything is good if the new MSA column is added, the renderer for the column is set, and the sequences within a cluster are of the same length.
---
{
  "order": 7
}