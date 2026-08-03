# Heat map tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Heat map

## Heatmap colors and global color scaling

1. Open properties, verify Heatmap Colors is checked
2. Uncheck Heatmap Colors
3. Check Heatmap Colors
4. Check Global Color Scaling
5. Uncheck Global Color Scaling

## Column label orientation

1. Open properties, verify Col Labels Orientation is set to Auto
2. Set Col Labels Orientation to Vert
3. Set Col Labels Orientation to Horz
4. Set Col Labels Orientation back to Auto

## Max heatmap columns

1. Open properties, verify Max Heatmap Columns is 100
2. Set Max Heatmap Columns to 3
3. Set Max Heatmap Columns to 1000
4. Set Max Heatmap Columns back to 100

## Show heatmap scrollbars

1. Open properties, verify Show Heatmap Scrollbars is checked
2. Uncheck Show Heatmap Scrollbars
3. Check Show Heatmap Scrollbars

## Is Heatmap toggle

1. Open properties, verify Is Grid is unchecked (viewer is in heatmap mode)
2. Check Is Grid — viewer switches to grid mode with row numbers and gridlines
3. Uncheck Is Grid — viewer returns to heatmap mode

## Filtering interaction

1. Open the filter panel
2. Set the AGE filter to range 20–40
3. Verify the heat map updates to show only filtered rows
4. Reset the AGE filter
5. Verify the heat map shows all rows again

## Table switching

> Note: requires spgi-100 dataset (open twice to get two tables for switching).
> Setup: Close all, open spgi-100.csv twice, go to the first table view, add Heat map, open Context Panel.

1. Set the Table property to the second spgi-100 table — heat map re-renders against it
2. Set the Table property back to the first spgi-100 table

## Range slider navigation

1. Drag the horizontal range slider to the left — content scrolls left
2. Drag the horizontal range slider to the right — content scrolls right
3. Double-click the horizontal range slider — resets to full range
4. Drag the vertical range slider up — content scrolls up
5. Drag the vertical range slider down — content scrolls down
6. Double-click the vertical range slider — resets to full range

## Column sorting

1. Double-click the AGE column header — rows sort ascending by AGE
2. Double-click the AGE column header again — rows sort descending
3. Double-click the AGE column header once more — sorting resets

## Selection interaction

1. Click a row in the heat map — it becomes the current row
2. Shift+drag to select a block of cells — selected rows are highlighted
3. Press Esc — selection clears

## Color scheme customization

1. Open properties, locate Linear Color Scheme
2. Change the linear color scheme to a different option
3. Open properties, locate Categorical Color Scheme
4. Change the categorical color scheme to a different option

## Draw every row

1. Open properties, verify Draw Every Row is unchecked
2. Check Draw Every Row
3. Scroll through the heat map
4. Uncheck Draw Every Row

## Layout save and restore

1. Set Col Labels Orientation to Vert
2. Check Global Color Scaling
3. Save the layout
4. Reset Col Labels Orientation to Auto and uncheck Global Color Scaling
5. Apply the saved layout — orientation and scaling restore
6. Delete the saved layout

## Layout saving with Is Heatmap toggle

> Note: uses spgi-100 dataset.
> Setup: open spgi-100.csv, add Heat map.

1. In the Property Pane, set Max Heatmap Columns to 200
2. Uncheck Is Heatmap — viewer switches to grid mode
3. Save the layout
4. Reset Max Heatmap Columns to 100 and check Is Heatmap
5. Apply the saved layout
6. Check Is Heatmap — re-enable heatmap mode
7. Delete the saved layout

---
{
  "order": 14,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
