1. Open SPGI.csv.
1. Click the `core` column header
1. Go to Context Pane > Chemistry > Rendering
1. Set Filter Type to `Categorical`
1. Using Filter icon open Filter Panel.
1. Check the Filter Panel:
1. The Structure tab should use Sketcher
1. The Primary Scaffold tab should contain categories (molecules).
1. Draw the c1cc2ccccc2cc1 molecule in the Structure tab Sketcher. Check the result.
1. Check the filtering by the categories from the Primary Scaffold tab.

***

1. Open the filter panel, clone view layout.
2. Set any filter on any view.
3. Save the layout.
3. Add some more filtering.
3. Apply the previously saved layout.

Expected: both layout and filtering of the dataframe should be applied

***

1. Open the Filter Panel.
2. Disable one filter.
3. Close the filter panel.
4. Open the filter panel - all other filters should be displayed on the panel.

***

1. Open the filter panel (FP).
1. View < Layout >clone view - check that the cloned Fp contains all filters.
2. Close any filter on the cloned FP.
3. Close the FP.
3. Open the FP - the closed filter should not be on the FP.
   
***

1. Add filtering by molecules, scaffold tree, scaffold tree viewer, categorical, numerical, scatterplot (filter by zoom), Bar chart and Pie chart (on click = filter), Pivot table (Row source = All)
2. Check the question mark - all the filtering should be listed
2. Click Reset Filter icon and check the question mark again - it should be empty, all the filtering should be reset

***

Reorder filter
1. Open SPGI.csv
1. Open the Filter Panel
1. In the hamburger menu of the Filter Panel, select Reorder filter.
1. A dialog opens. Check functionality 

***

1. Remove Structure filter.
2. Current value > Use as filter - check:
   * the molecule is sketched corectly,
   * the Structure filter is the first on the Filter Panel.

***

* [#1729](https://github.com/datagrok-ai/public/issues/1729)
* [#1747](https://github.com/datagrok-ai/public/issues/1747)
* [#1984](https://github.com/datagrok-ai/public/issues/1984)
* [#2165](https://github.com/datagrok-ai/public/issues/2165)
* [#2344](https://github.com/datagrok-ai/public/issues/2344)
* [#2358](https://github.com/datagrok-ai/public/issues/2358)
* [#1460](https://github.com/datagrok-ai/public/issues/1460)
* [#2185](https://github.com/datagrok-ai/public/issues/2185) Fuzzy filter enhancements
---
{
  "order": 16
}
