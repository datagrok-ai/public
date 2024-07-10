
## Default tooltip

### Grid: include visible columns in tooltip

1. Open a table with a few columns (all column names should be visible, e.g.,
   `DemoFiles/energy_uk.csv`)
2. Open grid properties and find `Show Visible Columns In Tooltip`
3. If it is unchecked (default), there should be no tooltip when hovering over grid cells
4. If it is unchecked (default), the tooltip should appear as you extend a column's width to push the last column(s) out of sight or extend the property panel to hide the last column(s).
5. Enable `Show Visible Columns In Tooltip`
6. Check that the tooltip is visible and remains the same both when all columns are visible and when some fall off the grid

### Viewers: uniform default tooltip

7. Open a scatter plot and a box plot
8. Hover over the viewers (including the main grid): the tooltip should consist of the same set of columns in the same order

### Viewers: actions in the context menu

9. Close all. Open demog dataset.
10. Open Histogram, Line Chart, Bar chart and Trellis plot.
11. Right-click on each viewer (including the main grid). The context menu under
   the `Tooltip` section should include:
   - `Hide`
   - `Edit...`
   - `Use as Group Tooltip`
   - `Remove Group Tooltip`

### Viewers: default tooltip visibility

12. Close all. Open demog dataset.
13. Open grid properties and enable `Show Visible Columns In Tooltip`
14. Open Scatter plot, Box plot, Histogram, Line Chart, Bar chart and Trellis plot.
15. Right-click on a viewer and select `Tooltip > Hide`
   - there should be no tooltip on hover over grid, scatter plot, box plot viewers
   - the tooltip for other viewers should not be affected
16. Right-click on a viewer and select `Tooltip > Show` (the option should appear instead of `Tooltip > Hide`)
   - the tooltip should be displayed on hover over all opened viewers


### Viewers: Edit tooltip

17. Close all. Open SPGI dataset.
18. Open grid properties and enable `Show Visible Columns In Tooltip`
19. Open a scatter plot and a box plot
20. Right-click on a viewer and select `Tooltip > Edit...`
   - dialog 'Set tooltip' should open
   - 'Search' checkbox should be checked, all variants should be picked (the current state: all columns are selected)
   - below the two checkboxes, there should be a grid listing the column type, name, and selection (three columns total)
   - try searching a column by its name in the search input over the grid
   - search is case-insensitive
   - under the grid, there are actions 'Reset group tooltip' and 'Design custom tooltip...'
   - the dialog footer includes the history icon, the buttons 'CANCEL' and 'OK'
21. Pick a few columns in the grid and click 'OK'
22. Hover over the viewers (including the main grid): the tooltip should consist of the same set of selected columns in the same order.


### Viewers: tooltip properties and heuristics for column selection

23. Open a table view
24. Add a scatter plot and a box plot
25. Find 'Show Tooltip' and 'Row Tooltip' in the 'Tooltip' section of the viewer properties
26. 'Row Tooltip' should be grayed-out (disabled), with an empty value
27. 'Show Tooltip' should have three options:
   - `inherit from table` (default)
   - `show custom tooltip`
   - `do not show`
28. Set 'Show Tooltip' to `show custom tooltip` for every viewer (including the grid, but its property 'Show Visible Columns In Tooltip' must be enabled)
29. Compare the custom tooltips with default table tooltip (add one more viewer for reference, as the main grid itself should be tested)
   - Without additional configuration, the custom viewer tooltip should use the same set of columns as the default tooltip

30. Right-click on a viewer with a custom tooltip and select `Tooltip > Hide`
   - The tooltip should be hidden only for this one viewer, the other viewers should use the tooltip as before, regardless of its type (custom or  default), e.g., if you have two scatter plots open and you hide a custom tooltip for one of them, it will disappear only for this one scatter plot.
31. Return the tooltip by choosing `Tooltip > Show` in the context menu
   - The tooltip should be visible and contain the same content as before
32. Check that an alternative way to toggle the custom tooltip visibility works:
   open viewer properties and switch 'Show Tooltip' from 'show custom tooltip' to 'do not show'
33. Hide the default tooltip (on a reference viewer from the previous steps)
   - The default tooltip should be hidden for all viewers that use it, while all custom viewer tooltips remain visible without any changes

