### PC plot

#### 1. Opening

1. Open SPGI, SPGI-linked1, SPGI-linked2
2. Go to **Tables > SPGI**
3. On the **Menu Ribbon**, click the **Add viewer** icon, and select **PC plot**
3. Close the viewer
1. Go to **Toolbox > Viewers > PC plot**

#### 2. Viewer basics 

1. Undock the viewer and move it around the screen
1. Dock the viewer in a different location
1. Open\Close the Context Panel
1. Display the help window
1. Close the viewer and return it by **Edit > Undo** (```CTRL + Z```)
1. Add a legend and change its position
1. Resize the viewer - **check mini-legend appearing and Self-adjustable viewer layout that keep viewers usable even in a small window**

#### 3. Context Panel

1. Open the **Context Panel**
5. Go to **Data > Table** and check tables switching (SPGI-linked2, SPGI-linked1, SPGI)
6. Go to **Value > Column Names** and change the list of selected columns - **Verify that axes are updated automatically**
1. Go to **Value > Log Columns** and set some columns to the log scale - **check the plot**
7. Go to **Data > Transformation**, set: 
   * Group by to `Whole Blood Assay Date 1`
   * Aggregate to `count(Id)`
   * Pivot to `Chemist 521`    
1. Save the Layout. Check

#### 4. Context menu

1. Go to **Context menu > Filter** and check/uncheck **Show Filters** - range slider on the viewer should appear/disappear
1. On the viewer, use range sliders and filter out some rows
1. Go to **Context menu > Filter** and check/uncheck **Show Filtered Out Lines** - verify that filtered out lines appear/disappear on the viewer
1.  Go to **Context menu > Y axis** and switch between Normalized and Global scales
1. Go to **Context menu > To Script** - a balloon with the script should appear
1. Go to **Context menu > General** and check all items
1. Go to **Context menu > Tooltip** and check all items
1. Open the Context Panel
1. Go to **Context menu > Columns / Selection / Properties** and verify that changes are consistent between the Context Panel and the context menu

#### 5. Color coding and the legend

1. Go to the **Color** info panel
2. Set Color to Chemical Space X
1. Go to **Context Panel > Style > Linear Color Scheme** and change it - **verify that scheme on the viewer is changed**
1. Grid > Chemical Space X: turn on the linear color coding and change the color scheme - **verify that scheme on the viewer is changed**
1. Grid > Chemical Space X: switch to the conditional color coding - **verify that legend appears and scheme on the viewer is changed**
1. Go to **Context Panel > Legend** and check Legend visibility and position settings
3. Save to Layout. Check

#### 6. Pick Up / Apply  

1. Add the second PC plots
1. For the first plot:
  * Change the set of axes
  * Switch some axes to the log scale
  * Enable the categorical legend and change its position
  * On the legend, change some colors
  * Add the title
1. Right-click the first PC plot. Select **Pick up/Apply > Pick up** 
1. Right-click the second PC plot. Select **Pick up/Apply > Apply**
1. Change the axes on the first PC plot - **the second plot should not be affected**
1. Adjust the range slider on the second PC plot - **the first plot should update to show the filtered lines, but its own range sliders should remain unchanged and stay in their original positions**

#### 7. Layout and Project save/restore
1. Save the layout
2. Add some more viewers
1. Apply the saved layout - **verify that only two PC plots are displayed**
1. Save the project
1. Close All
1. Open the saved project - **check**

#### 8. Filtering

1. Close All
1. Open SPGI
1. Add a PC plot
1. Use range sliders on the axes to filter data
1. Open the **Filter Panel** and apply some filters - **Filtering should respect both Filter Panel and PC plot filters**
1. Go to PC plot, **Context menu > Reset view** (or double-click the white space in the viewer) - **Only the filtering on the PC plot should be reset, and the range sliders should return to their default positions**
1. Use range sliders on the axes to filter data one more time
1. On the **Filter Panel** click **Reset filter** - this should reset all applied filters and restore the range sliders on the PC plot to their default positions

---
{
  "order": 8,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}