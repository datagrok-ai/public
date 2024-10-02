### Viewers: tooltip properties and heuristics for column selection

1. Open a table view
2. Add a scatter plot and a box plot
3. Find 'Show Tooltip' and 'Row Tooltip' in the 'Tooltip' section of the viewer properties
4. 'Row Tooltip' should be grayed-out (disabled), with an empty value
5. 'Show Tooltip' should have three options:
    - `inherit from table` (default)
    - `show custom tooltip`
    - `do not show`
6. Set 'Show Tooltip' to `show custom tooltip` for every viewer (including the grid, but its property 'Show Visible
   Columns In Tooltip' must be enabled)
7. Compare the custom tooltips with default table tooltip (add one more viewer for reference, as the main grid itself
   should be tested)
    - Without additional configuration, the custom viewer tooltip should use the same set of columns as the default
      tooltip

8. Right-click on a viewer with a custom tooltip and select `Tooltip > Hide`
    - The tooltip should be hidden only for this one viewer, the other viewers should use the tooltip as before,
      regardless of its type (custom or default), e.g., if you have two scatter plots open and you hide a custom tooltip
      for one of them, it will disappear only for this one scatter plot.
9. Return the tooltip by choosing `Tooltip > Show` in the context menu
    - The tooltip should be visible and contain the same content as before
10. Check that an alternative way to toggle the custom tooltip visibility works:
    open viewer properties and switch 'Show Tooltip' from 'show custom tooltip' to 'do not show'
11. Hide the default tooltip (on a reference viewer from the previous steps)
    - The default tooltip should be hidden for all viewers that use it, while all custom viewer tooltips remain visible
      without any changes

---
{
"order": 6
}
