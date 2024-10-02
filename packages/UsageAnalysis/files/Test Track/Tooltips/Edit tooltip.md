### Viewers: Edit tooltip

1. Close all. Open SPGI dataset.
2. Open grid properties and enable `Show Visible Columns In Tooltip`
3. Open a scatter plot and a box plot
4. Right-click on a viewer and select `Tooltip > Edit...`
   - dialog 'Set tooltip' should open
   - 'Search' checkbox should be checked, all variants should be picked (the current state: all columns are selected)
   - below the two checkboxes, there should be a grid listing the column type, name, and selection (three columns total)
   - try searching a column by its name in the search input over the grid
   - search is case-insensitive
   - under the grid, there are actions 'Reset group tooltip' and 'Design custom tooltip...'
   - the dialog footer includes the history icon, the buttons 'CANCEL' and 'OK'
5. Pick a few columns in the grid and click 'OK'
6. Hover over the viewers (including the main grid): the tooltip should consist of the same set of selected columns in the same order.


---
{
  "order": 5
}