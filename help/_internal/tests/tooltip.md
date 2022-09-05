# Tests: tooltip

## Default tooltip

### Grid: include visible columns in tooltip

1. Open a table with a few columns (all column names should be visible, e.g.,
   `DemoFiles/energy_uk.csv`)
1. Open grid properties and find `Show Visible Columns In Tooltip`
1. If it is unchecked (default), there should be no tooltip when hovering over
   grid cells
1. If it is unchecked (default), the tooltip should appear as you extend a
   column's width to push the last column(s) out of sight or extend the property
   panel to hide the last column(s). The tooltip should include only columns
   that cannot be seen without scrolling to the right
1. Enable `Show Visible Columns In Tooltip`
1. Check that the tooltip is visible and remains the same both when all columns
   are visible and when some fall off the grid

### Viewers: uniform default tooltip

1. Open a table view
1. Open grid properties and enable `Show Visible Columns In Tooltip`
1. Open a scatter plot and a box plot (for now, skip PC plot)
1. Hover over the viewers (including the main grid): the tooltip should consist
   of the same set of columns in the same order

### Viewers: actions in the context menu

1. Open a table view
1. Open a viewer of each type
1. Right-click on each viewer (including the main grid). The context menu under
   the `Tooltip` section should include:
   - `Hide`
   - `Set Default Tooltip`
   - `Set ${Viewer} Tooltip` (only for grid, scatter plot, box plot; may be
     added for PC plot)
   - `Use as Group Tooltip`
   - `Remove Group Tooltip`

### Viewers: default tooltip visibility

1. Open a table view
1. Open grid properties and enable `Show Visible Columns In Tooltip`
1. Open a viewer of each type
1. Right-click on a viewer and select `Tooltip > Hide`
   - there should be no tooltip on hover over
     - grid
     - scatter plot
     - box plot
     - heat map
   - the tooltip for other viewers should not be affected
1. Right-click on a viewer and select `Tooltip > Show` (the option should appear
   instead of `Tooltip > Hide`)
   - the tooltip should be displayed on hover over
     - grid
     - scatter plot
     - box plot
     - heat map
   - the tooltip for other viewers should not be affected

### Viewers: adjust columns in default tooltip

1. Open a table view
1. Open grid properties and enable `Show Visible Columns In Tooltip`
1. Open a scatter plot and a box plot (for now, skip PC plot)
1. Right-click on a viewer and select `Tooltip > Set Default Tooltip`
   - dialog 'Set default tooltip' should open
   - 'Show tooltip' checkbox should be checked
   - 'All on' checkbox should be checked (*or should it be changed to reflect
    the current state: all columns are selected in the grid below or not?*)
   - below the two checkboxes, there should be a grid listing the column type,
     name, and selection (three columns total)
   - try searching a column by its name in the search input over the grid
     - search is case-insensitive
     - the results include not only the exact matches (substring search works)
   - under the grid, there are actions 'Reset group tooltip' and 'Sketch
     form...'
   - the dialog footer includes the history icon, the buttons 'CANCEL' and 'OK'
1. Pick a few columns in the grid and click 'OK'
1. Hover over the viewers (including the main grid): the tooltip should consist
   of the same set of selected columns in the same order

### Viewers: form tooltip

1. Open a table view
1. Open a scatter plot and a box plot (for now, skip PC plot)
1. Right-click on a viewer and select `Tooltip > Set Default Tooltip`
1. Select `Sketch form...` in the 'Set default tooltip' dialog
1. Make changes to a displayed form and click 'Close and apply'
1. Hover over the viewers (including the main grid): the tooltip should
   correspond to the template sketched in the previous step

### Viewers: switch from form to columns in default tooltip

Preconditions: follow the steps in 'Viewers: form tooltip'

1. Right-click on a viewer and select `Tooltip > Set Default Tooltip`
1. Select tooltip columns in the dialog and click 'OK'
   - the default tooltip is updated so that selected columns are shown
   - the tooltip is the same for all viewers

## Viewer-specific tooltip

### Scatter plot: set tooltip columns

1. Open a table view
1. Add a scatter plot
1. Right-click on the scatter plot and select `Tooltip > Set Scatter plot
   Tooltip`
1. Set 'Show tooltip' to 'show custom tooltip'
1. Select tooltip columns in the dialog and click 'OK'
   - Scatter plot tooltip is updated so that selected columns are shown
   - The main grid uses the default tooltip (as well as other scatter plot
     instances, open a new one to check)

### Scatter plot: custom form tooltip

1. Open a table view
1. Add a scatter plot
1. Right-click on the scatter plot and select `Tooltip > Set Scatter plot
   Tooltip`
1. Set 'Show tooltip' to 'show custom tooltip'
1. Select `Sketch form...` in the 'Set custom tooltip' dialog
1. Make changes to a displayed form and click 'Close and apply'
1. Hover over the scatter plot
   - The tooltip should correspond to the template sketched in the previous step
   - Other scatter plots should use the default tooltip (open a new plot to
     check)

### Scatter plot: switch from form default tooltip to custom tooltip with columns

1. Follow the steps in 'Viewers: form tooltip'
1. Follow the steps in 'Scatter plot: set tooltip columns' in the same table
   view
   - Scatter plot tooltip is updated so that selected columns are shown
   - Default tooltip shows the form

### Scatter plot: switch from form default tooltip to custom form tooltip

1. Follow the steps in 'Viewers: form tooltip'
1. Follow the steps in 'Scatter plot: custom form tooltip' in the same table
   view
   - In a sketch view for scatter plot tooltip, the provided template should be
   default (the state is refreshed and doesn't copy the form in the default
   tooltip)
   - Scatter plot tooltip is updated so that the form is shown
   - Default tooltip shows the old form (not the form configured for the scatter
     plot)

### Viewers: tooltip properties and heuristics for column selection

1. Open a table view
1. Add a scatter plot and a box plot (for now, skip PC plot)
1. Find 'Show Tooltip' and 'Row Tooltip' in the 'Tooltip' section of the viewer
   properties (including the grid)
1. 'Row Tooltip' should be grayed-out (disabled), with an empty value
1. 'Show Tooltip' should have three options:
   - `inherit from table` (default)
   - `show custom tooltip`
   - `do not show`
1. Set 'Show Tooltip' to `show custom tooltip` for every viewer (including the
   grid, but its property 'Show Visible Columns In Tooltip' must be enabled)
1. Compare the custom tooltips with default table tooltip (add one more viewer
   for reference, as the main grid itself should be tested)
   - Without additional configuration, the custom viewer tooltip should use the
     same set of columns as the default tooltip

### Viewers: custom tooltip visibility

1. Follow the steps in 'Viewers: tooltip properties and heuristics for column
   selection'
1. Right-click on a viewer with a custom tooltip and select `Tooltip > Hide`
   - The tooltip should be hidden only for this one viewer, the other viewers
     should use the tooltip as before, regardless of its type (custom or
     default), e.g., if you have two scatter plots open and you hide a custom
     tooltip for one of them, it will disappear only for this one scatter plot
1. Return the tooltip by choosing `Tooltip > Show` in the context menu
   - The tooltip should be visible and contain the same content as before
1. Check that an alternative way to toggle the custom tooltip visibility works:
   open viewer properties and switch 'Show Tooltip' from 'show custom tooltip'
   to 'do not show'
1. Hide the default tooltip (on a reference viewer from the previous steps)
   - The default tooltip should be hidden for all viewers that use it, while all
     custom viewer tooltips remain visible without any changes

### Scatter plot: switch from custom tooltip to default tooltip with columns

1. Follow the steps in 'Scatter plot: set tooltip columns' or 'Scatter plot:
   custom form tooltip'
1. Follow some of the steps in 'Viewers: adjust columns in default tooltip' to
   switch to the default tooltip from the scatter plot context menu
   - The scatter plot now uses the default tooltip, which is updated so that
     the selected columns are shown
   - Default tooltip for other viewers looks the same as for the scatter plot

### Scatter plot: switch from custom tooltip to default form tooltip

1. Follow the steps in 'Scatter plot: set tooltip columns' or 'Scatter plot:
   custom form tooltip'
1. Follow some of the steps in 'Viewers: form tooltip' to switch to the default
   tooltip from the scatter plot context menu
   - In a sketch view, the provided template should be default
   - The scatter plot now uses the default tooltip, which is updated so that
     the form is shown
   - Default tooltip for other viewers looks the same as for the scatter plot
