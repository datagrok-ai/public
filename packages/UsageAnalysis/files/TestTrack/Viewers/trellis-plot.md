# Trellis plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Trellis plot

## Inner viewer types

1. Switch inner viewer to **Scatter plot** using the viewer type selector at the top
2. Set X to WEIGHT, Y to HEIGHT, Color to RACE
3. Switch to **Bar chart**
4. Set Split to RACE, Value to AGE, Value Aggr Type to avg
5. Switch to **Histogram**
6. Set Value to AGE, Split to RACE
7. Switch to **Line chart**
8. Set X to STARTED, Y to AGE, Split to RACE
9. Switch to **Box plot**
10. Set Category to SEX, Value to AGE
11. Switch to **Pie chart**
12. Set Category to RACE
13. Switch to **Density plot**
14. Set X to WEIGHT, Y to HEIGHT
15. Switch to **Summary**
16. Set Visualization to bars
17. Switch to **Sparklines**
18. Set Sparkline Type to Bar Chart
19. Switch to **PC plot**
20. Set Color to SEX

## Global scale

1. Switch inner viewer to **Scatter plot**
2. Open Context Panel > **Axes > Global Scale** -- enable it
3. Disable **Global Scale**
4. Re-enable **Global Scale**

## Axes visibility

1. Switch inner viewer to **Scatter plot**
2. Open Context Panel > **Axes > Show X Axes** > set to **Always**
3. Set to **Never**
4. Set to **Auto**
5. Repeat for **Show Y Axes**: Always, Never, Auto
6. Toggle **Show Range Sliders** off
7. Toggle **Show Range Sliders** on

## Range sliders with global scale

1. Switch inner viewer to **Scatter plot**
2. Enable **Global Scale** and **Show Range Sliders**
3. Set **Show X Axes** to **Always**
4. Hover over the X axis area -- a range slider should appear
5. Drag the range slider to narrow the visible data range -- all inner viewers should update
6. Right-click on the trellis plot and select **Reset Inner Range Sliders** -- sliders reset to full range
7. Repeat for Y axis (set **Show Y Axes** to **Always**)

## Gridlines

1. Open Context Panel > **Show Gridlines** > set to **always**
2. Set to **never**
3. Set to **auto**

## Tiles mode

1. Open Context Panel > **Tiles** -- enable tiles view
2. Set **Tiles Per Row** to 2
3. Set **Tiles Per Row** to 6
4. Disable **Tiles**

## Category management

1. Set X to SEX, Y to RACE
2. Add a second X column DIS_POP
3. Remove the second X column
4. Toggle **Show X Labels** off
5. Toggle **Show Y Labels** off
6. Re-enable both

## Pack categories

1. Set X to SEX, Y to RACE
2. Open the filter panel and filter out one RACE category (e.g., uncheck "Asian")
3. Verify **Pack Categories** is enabled -- the empty row should disappear
4. Disable **Pack Categories** -- the empty row reappears
5. Re-enable **Pack Categories**

## On Click functionality

1. Open Context Panel > **On Click** > set to **Select**
2. Click any non-empty trellis cell -- rows matching that cell's categories should become selected
3. Change the inner viewer type -- selection should NOT change
4. Click another non-empty cell -- selection should change
5. Change any axis -- selection should NOT change
6. Set **On Click** to **Filter**
7. Click any non-empty cell -- only matching rows should remain visible
8. Change the inner viewer -- filtering should NOT change
9. Click another non-empty cell -- filtering should update
10. Change any axis -- filtering should be reset
11. Click a non-empty cell
12. Press **ESC** -- filtering and selection should reset
13. Click a non-empty cell -- some rows should be filtered
14. Open the **Filter Panel** and apply additional filters -- filtering should respect both filter panel and trellis filters
15. Right-click > **Reset view** -- only the trellis filtering should reset
16. Set **On Click** to **None** -- clicking any cell should not change filtering or selection

## Selectors

1. Toggle **Show X Selectors** off
2. Toggle **Show Y Selectors** off
3. Toggle **Show Control Panel** off
4. Re-enable all three

## Allow viewer full screen

1. Hover over a trellis cell -- a full-screen icon should appear
2. Click the full-screen icon -- the inner viewer should expand
3. Close the full-screen view
4. Disable **Allow Viewer Full Screen** in Context Panel -- full-screen icon should no longer appear on hover

## Scrolling

1. Set X to SITE -- many categories should appear
2. Verify horizontal scroll slider appears at the bottom
3. Drag the horizontal slider to scroll through categories
4. Set Y to RACE and X to SITE -- verify vertical scroll also appears
5. Use mouse wheel over the trellis to scroll vertically

## Legend

1. Switch inner viewer to **Scatter plot** and set a color column to SEX
2. Set **Legend Visibility** to **Always**
3. Set **Legend Position** to Left, Right, Top, Bottom
4. Set **Legend Visibility** to **Never**

## Context menu

1. Right-click on a trellis cell -- context menu should include a group named after the inner viewer type (e.g., **Pie chart**)
2. The inner viewer group should contain the inner viewer's context menu items
3. Standard viewer menu items (Properties, Clone, etc.) should appear below

## Inner viewer properties

1. Switch inner viewer to **Scatter plot**
2. Click the gear icon on the trellis plot title bar to open properties
3. Switch to the inner viewer tab at the top of the Context Panel
4. Change an inner viewer property (e.g., change X or Y column)
5. Verify the change applies to all cells

## Use in Trellis

1. Close the trellis plot
2. Add a **Scatter plot**, configure it (set X to AGE, Y to HEIGHT, color to SEX)
3. Right-click on the scatter plot > **General** > **Use in Trellis** -- a trellis plot should appear with the scatter plot as inner viewer preserving settings
4. Close the trellis plot and the scatter plot
5. Add a **Bar chart**, configure it
6. Right-click > **General** > **Use in Trellis** -- verify trellis appears with bar chart
7. Close the trellis plot and the bar chart
8. Add a **Histogram**, configure it
9. Right-click > **General** > **Use in Trellis** -- verify
10. Close the trellis plot and the histogram
11. Add a **Line chart**, configure it
12. Right-click > **General** > **Use in Trellis** -- verify
13. Close the trellis plot and the line chart
14. Add a **Box plot**, configure it
15. Right-click > **General** > **Use in Trellis** -- verify

## Auto layout

1. Enable **Auto Layout** (should be enabled by default)
2. Resize the trellis plot to a small size -- control panel, labels, and selectors should hide automatically
3. Resize back to a large size -- controls reappear
4. Disable **Auto Layout** -- all controls remain visible regardless of size

## Title and description

1. Open Context Panel > **Description > Show Title** -- enable it
2. Set **Title** to "My Trellis"
3. Set **Description** to "Test description"
4. Change **Description Position** to Bottom, Top, Left, Right

## Label orientation

1. Set X to RACE
2. Open Context Panel > **X Labels Orientation** -- change to Horizontal
3. Change to Vertical
4. Change to Auto
5. Repeat for **Y Labels Orientation**

## Pick Up / Apply

3. Add two trellis plots
4. For the first trellis plot: set Y axis to R1, switch inner viewer to bar chart, enable legend and change its position, add a title
5. Right-click the first trellis plot > **Pick up/Apply > Pick up**
6. Right-click the second trellis plot > **Pick up/Apply > Apply** -- second plot should match the first
7. Change the X or Y axis on the first trellis plot -- the second plot should not be affected
8. Adjust the range slider on the second trellis plot -- the first plot should not be affected

## Layout and Project save/restore

1. Save the layout
2. Add some more viewers
3. Apply the saved layout -- verify that only original viewers are displayed
4. Save the project
5. Close all
6. Open the saved project -- verify the correct viewers are displayed


## Viewer filter formula

1. Open Context Panel > **Data > Filter**
2. Set filter formula to `${AGE} > 40`
3. Clear filter formula

## Multi Curve inner viewer (and table switching)

1. Open dataset from **Files > Demo > curves.csv**
2. Go back to the demog view
3. For trellis plot, go to the Context Panel > Data  and set Table to curves.
4. On the trellis plot, select **Multi Curve viewer** as inner viewer
5. Set X and Y axes on the viewer
6. Change number of categories on the axes using +/-
7. Move/resize zoom slider for X and Y axes
8. Click the **Gear** icon and check properties

## To Script

1. Right-click the trellis plot and select To Script
2. Verify a balloon with the script appears
3. Close the balloon
## Keyboard navigation

1. Click on a trellis cell to select it
2. Press **Right Arrow** -- selection should move to the next cell
3. Press **Left Arrow** -- selection moves back
4. Press **Down Arrow** -- selection moves to the row below
5. Press **Up Arrow** -- selection moves to the row above
6. Press **ESC** -- trellis filter should reset
---
{
  "order": 6,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/curves.csv"]
}
