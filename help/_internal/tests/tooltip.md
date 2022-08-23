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
1. Open a scatter plot and a box plot (for now, skip PC plot and heat map)
1. Hover over the viewers (including the main grid): the tooltip should consist
   of the same set of columns in the same order

### Viewers: actions in the context menu

1. Open a table view
1. Open a viewer of each type
1. Right-click on each viewer (including the main grid). The context menu under
   the `Tooltip` section should include:
   - `Hide`
   - `Set Default Tooltip`
   - `Set ${Viewer} Tooltip` (only for scatter plot, box plot; may be added for
     PC plot, heat map))
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
1. Open a scatter plot and a box plot (for now, skip PC plot and heat map)
1. Right-click on a viewer and select `Tooltip > Set Default Tooltip`
   - dialog 'Select tooltip columns' should open
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
1. Open a scatter plot and a box plot (for now, skip PC plot and heat map)
1. Right-click on a viewer and select `Tooltip > Set Default Tooltip`
1. Select `Sketch form...` in the 'Select tooltip columns' dialog
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
1. Select tooltip columns in the dialog and click 'OK'
   - Scatter plot tooltip is updated so that selected columns are shown

### Scatter plot: custom form tooltip

1. Open a table view
1. Add a scatter plot
1. Right-click on the scatter plot and select `Tooltip > Set Scatter plot
   Tooltip`
1. Select `Sketch form...` in the 'Select tooltip columns' dialog
1. Make changes to a displayed form and click 'Close and apply'
1. Hover over the scatter plot
   - the tooltip should correspond to the template sketched in the previous step
   - other scatter plots should use the default tooltip (open a new plot to
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
   - Scatter plot tooltip is updated so that the form is shown
   - Default tooltip shows the old form (not the form configured for the scatter
     plot)
