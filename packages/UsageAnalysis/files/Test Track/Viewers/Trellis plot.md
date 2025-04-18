### Trellis plot

#### 1. Opening
1. Open demog.csv
2. On the **Menu Ribbon**, click the **Add viewer** icon, and select **Trellis plot**
3. Close the viewer
1. Go to **Toolbox > Viewers > Trellis plot**
---

#### 2.  Viewer basics 

1. Undock the viewer and move it around the screen
1. Dock the viewer in a different location
1. Switch between inner viewers:
   * Change their settings on the top of the trellis plot
   * Hover and select elements inside the viewer
1. Open\Close the Context Panel
1. Display the  help window
1. Close the viewer and return it by **Edit > Undo** (```CTRL + Z```)
---
   
#### 3. Context Panel

1. Open the **Context Panel** - check its appearance and UI
1. Undock the panel, move it around the workspace and dock again
1. Check search on the Context Panel
1. Change some settings values - **the viewer should be updated accordingly**
1. Add or remove the X (Y) column directly from the trellis plot
1. Verify the number of X (Y) columns in the Context Panel - **the number of columns in the Context Panel should be immediately updated when a column is added or removed**
---

#### 4. Context menu

1. Switch the inner viewer and check its specific item in the context menu
1. Go to **Context menu > To Script** - a balloon with the script should appear
1. Go to **Context menu > General** and check all items
1. Go to **Context menu > Tooltip** and check all items
1. Open the Context Panel
1. Go to **Context menu > Properties** and verify that property changes are consistent between the Context Panel and the context menu
---

#### 5. OnClick functionality (check for each inner viewer)

1. Go to **Context menu > On Click > Select**
1. Click any empty cell - **there should be no selected rows**
1. Click any non-empty cell - **check selection**
1.  Change the inner viewer  - **selection should NOT change**
1. Click non empty cell - **selection should change**
1. Change any axis - **selection should NOT change**
1. Go to **Context menu > On Click > Filter** - all selected rows should be filtered, Rows Source should be changed to All
1. Click any empty cell - **all rows should be filtered out**
1. Click any non-empty cell - **check filtering**
1. Change the inner viewer - **filtering should NOT change**
1. Click non empty cell - **check filtering**
1. Change any axis - **filtering should be reset**
1. Click non empty cell
1. Press ESC - **this should reset filtering and selection**
1. Click non empty cell - **some rows should be filtered**
1. Open the **Filter Panel** and apply some filters - **Filtering should respect both Filter Panel and trellis plot filters**
1. Go to **Context menu > Reset view** - only the filtering of the trellis plot should be reset
1. Go to **Context menu > On Click > None** - clicking any cell should not change filtering or selection
---

#### 6. Pick Up / Apply  ([#2116](https://github.com/datagrok-ai/public/issues/2116))

1. Close All
1. Open SPGI 
1. Add two trellis plots
1. For the first trellis plot:
   * Set the Y axis to R1
   * Switch the inner viewer to bar chart
   * Enable the legend and change its position
   * Add the title
1. Right-click the first trellis plot. Select **Pick up/Apply > Pick up** 
1. Right-click the second trellis plot. Select **Pick up/Apply > Apply**
1. Change the X or Y axis on the first trellis plot - **the second plot should not be affected**
1. Adjust the range slider on the second trellis plot using +/- on the axes, scroll - **the first plot should not be affected**
---

#### 7. Layout and Project save/restore
1. Save the layout
2. Add some more viewers
1. Apply the saved layout - **verify that only two trellises are displayed**
1. Save the project
1. Close All
1. Open the saved project - **verify that only two trellises are displayed**
---

#### 8. Floating viewer after applying layout

1. Close All
1. Open demog.csv
1. Use `Ctrl + '-'` (or `Cmd + '-'` on Mac) to zoom out the screen view
1. Add a trellis plot 
1. Undock the viewer and drag it to the bottom of the screen
1. Save the layout.
1. Reset Screen Zoom: Zoom in to return the screen to its original size
1. Apply the saved layout

**Expected Results:**

  * The layout should automatically adjust to fit the current screen size.
  * All viewers should remain accessible and properly positioned on the screen.
---
{
  "order": 13,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
