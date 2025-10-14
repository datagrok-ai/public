### Filter panel

- Open the filter panel, clone view layout (Top Menu > View > Layout > Clone view).
- Set any filter on any view, it should apply filter to both views.
- Save the layout (on Toolbox: Layout > Save).
- Add some more filtering.
- Apply the previously saved layout.

Expected result: both layout and filtering of the dataframe should be applied to the dataset from layout

- Add filtering by molecules, categorical, numerical: Filter Panel > Hamburger menu > Add Filter > Scaffold Tree Filter
- Add viewers and apply filtering on them: scaffold tree viewer, scatterplot (filter by zoom), Bar chart, Pie chart (on click = filter), Pivot table (Row source = All)
- Check the question mark on Filter Panel - all the filtering should be listed
- Click Reset Filter icon (looks like reload) and check the question mark again - it should be empty, all the filtering should be reset
- No hidden columns should be visible in filter panel
---
{
"order": 1
}