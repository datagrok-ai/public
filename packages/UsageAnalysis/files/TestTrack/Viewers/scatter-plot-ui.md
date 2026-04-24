Manual checklist. Not included in Playwright automation.

### Tooltip

1. On the viewer, right-click the canvas and go to Tooltip — check AGE and WEIGHT in the column list
2. Right-click again and set Tooltip > Show Tooltip to Never
3. Go to the Context Panel > Tooltip, set Show Tooltip to Always
4. Save the layout, close the scatter plot, apply the layout — verify tooltip settings are preserved

### Title, description and range sliders

1. Go to the Context Panel > Misc, set Title to "Test Plot", set Description to "Test description"
2. Set Title Position to Top, Description Position to Bottom
3. Check the range sliders on the X and Y axes — drag them to narrow the visible range
4. Save the layout, close the scatter plot, apply the layout — verify title, description, and range slider positions are preserved

### Formula lines

1. On the viewer, right-click the canvas and select Tools > Formula Lines
2. In the Formula Lines editor, click ADD NEW and add a horizontal line, a vertical line, and a formula line
3. Switch between the Viewer and Dataframe tabs — verify Viewer lines appear only on this viewer
4. Save the layout, close the scatter plot, apply the layout — verify only Viewer-tab formula lines are restored
5. Delete all created lines

### Table switching and filter expression

1. Close all and open SPGI, SPGI-linked1, SPGI-linked2
2. Add a scatter plot
3. Go to the Context Panel > Data, switch Table to SPGI-linked1 — verify axes update
4. Switch Table to SPGI-linked2 — verify axes update
5. Switch Table back to SPGI
6. In the Context Panel > Data, set Filter to `${CAST Idea ID} <636500` — verify points are filtered on the plot

### Zoom and filter combinations with jitter

1. On the viewer, set X to AGE, Y to HEIGHT
2. Go to the Context Panel > Marker, set Jitter Size to 30
3. Zoom in on the scatter plot using mouse wheel
4. Go to the Context Panel > Data, set Zoom and Filter to 'zoom by filter'
5. Open the filter panel and filter RACE to 'Asian' — verify the plot zooms to fit filtered points
6. Set Zoom and Filter to 'pack and zoom by filter' — verify points repack and zoom adjusts
7. Drag the X range slider while jitter is active — verify points and filtering remain consistent
8. Set Jitter Size back to 0

### Grid color coding sync with scatter plot

1. Go to the grid, right-click the header of SEX and select Color coding > Categorical
2. Go to the scatter plot, set Color to SEX — verify colors match the grid's categorical scheme
3. Go to the grid, right-click the header of AGE and select Color coding > Linear, Edit — change the color scheme and invert it
4. On the scatter plot, set Color to AGE — verify the scatter plot uses the same inverted scheme
5. Go to the grid, right-click the header of AGE and select Color coding > Edit — change the scheme again
6. Verify the scatter plot updates to match

### Markers: structure and datetime columns, Markers Map

1. Close all and open SPGI
2. Add a scatter plot
3. Go to the Context Panel > Marker, set Markers to the Structure column — verify structure markers render on the plot
4. Set Markers to the Synthesis Date column (datetime) — verify datetime markers render
5. Change Markers Map values for the datetime column — verify the legend and plot remain consistent
6. Set Markers to Primary Series Name (categorical) — verify standard shape markers

### Selection with jitter changes and deselection

1. On the viewer, set X to AGE, Y to HEIGHT
2. Go to the Context Panel > Marker, set Jitter Size to 20, Jitter Size Y to 15
3. Hold Shift and drag a rectangle to select points — verify selection in the grid
4. Change Jitter Size to 30 — verify the selection is preserved in the grid
5. Reset Jitter Size Y to 0 while Jitter Size remains 30 — verify the selection is still preserved
6. Hold Ctrl+Shift and drag a rectangle over some selected points — verify those points are deselected
7. Press Escape to deselect all

### Structure as X axis

1. Close all and open SPGI
2. Add a scatter plot
3. On the viewer, click the X column selector and set X to Structure — verify structures render along the X axis
4. Save the layout, close the scatter plot, apply the layout — verify Structure axis is preserved