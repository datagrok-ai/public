Open linked datasets

1. Add a structure filter (draw or current value > use as filter).
2. Disable the structure filter.
3. Close the filter panel.
3. Open  the filter panel.
3. Turn on the structure filter.

***

1. Set a filter.
2. Close the filter panel.
3. Use Current value > Use as filter.

***

1. Set a filter.
2. Close the filter panel.
3. Go to the mol column's hamburger menu > filter.
4. Draw another structure.
5. Click Add filter.

***

1. Add filters.
2. View > Layout > Clone View.
3. Disable filters on the clone.
3. Close-open the filter panel.
4. Enable filters.
5. Clear filters on clone.

***

1. Clone a view.
2. Open the FP on the original view.
3. Open the FP on the clone.
2. Add a structure filter on the clone.
3. Close the filter on the clone.
3. Go to the original view and close the FP - the highlight shouldn't remain.

***

1. Remove Structure filter.
2. Current value > Use as filter - check:
   * the molecule is sketched corectly,
   * the Structure filter is the first on the Filter Panel.

***

Chem: Substructure Filter: Not terminated on filter disable:
1. Open tests/smi10K.csv dataset or other dataset with more than 10K structures
2. Open filter panel and start the structure search.
3. While search is ongoing, disable filter.

Expected result: search is terminated.

***

1. Open linked datasets.
2. Open filters and make sure some structure column is added to the filter panel.
3. Add another view with the same table (e.g. clone current view).
4. Open filters and make sure same structure column is added to the filter panel in the 2nd view.
4. Right-click some structure > Current value > Use as filter.

Expected result: filter is applied, filters are synchronized in two open views.

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
