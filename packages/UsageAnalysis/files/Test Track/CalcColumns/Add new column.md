1. Open:
   * Table form local storage
   * Table form Home dir
   * query result
   * GetTop100 or GetAll result
1. Click the **Add new column** icon.
1. Add a formula using columns from the dataset.
1. Add another new column using the column added in the previous step.
1. In the grid, change the used columns names and some values - newly added columns should change their formulas and values accordingly.
1. Save a project with datasync (where available).
1. Close All
1. Open the project (saved with datasync).
1. Change the name of the column used in formula - it should be changed in formula as well.
1. Change any values in the column used in formula - values in the resulting column should be updated accordingly.

Expected: calc column is added and its values change according to the changes in table columns.

---
{
  "order": 1
}
