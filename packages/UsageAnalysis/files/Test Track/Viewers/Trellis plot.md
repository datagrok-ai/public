### Trellis plot

1. Open the **demog** dataset.
2. On the **Menu Ribbon**, click the **Add viewer** icon, and select **Trellis plot**. 
   * **_Expected result:_** **Trellis plot** visualized viewer opens without errors in the console. 
3. Close previously opened viewer. On the **Viewers tab**, click **Trellis plot** icon. 
   * **_Expected result:_** **Trellis plot** visualized viewer opens without errors in console. 
4. Working with a viewer:
    * Moving the viewer by the workspace
    * Navigation inside the viewer
    * Hover and select elements inside the viewer
    * Open\Close the Context Panel
    * Display of help window
    * Return viewer after closing by **Edit | Undo** (```CTRL + Z```)
    * Interact with all viewer`s elements directly on the viewer.
   * **_Expected Result_**: Changes should be reflected on the viewer with no errors. 
5. On the **Trellis plot** viewer, click the **Gear** icon. The **Context Panel** opens. Check:
    * Appearance and UI. Context Panel`s style shouldn't be broken
    * Detach the panel and move it around the workspace
    * Search on the Context Panel
    * Change some settings values
   * **_Expected Result:_** Changed settings should be reflected in the Trellis plot viewer without errors, and the viewer should be updated accordingly.
6. **Context menu** (focus on options: inner viewer changes, properties, Clone, save as PNG, Save to Gallery, Embed, Full Screen, Close)
   * **_Expected Result:_** Changes should be reflected in the viewer without errors, and the viewer should be updated accordingly.
6. **On Click** functionality:
   * Right-click the viewer and select **On Click > Select**
   * Click any empty cell - there should be no selected rows
   * Click any non-empty cell - check selection
   * Change the inner viewer  - selection should NOT change
   * Click non empty cell - selection should change
   * Change any axis - selection should NOT change
   * Right-click the viewer and select **On Click > Filter** - all selected rows should be filtered, Rows Source should be changed to All
   * Click any empty cell - all rows should be filtered out
   * Click any non-empty cell - check filtering
   * Change the inner viewer - filtering should NOT change
   * Click non empty cell - check filtering
   * Change any axis - filtering should be reset
   * Click non empty cell
   * Press ESC - this should reset filtering and selection
   * Right-click the viewer and select **On Click > None** - clicking any cell should not change filtering of selection
7. **Tooltip** testing (main focus on options: Hide, Edit, Use a group tooltip, use as group tooltip, remove group tooltip)
   * Changes to the properties should be reflected in the Trellis plot viewer without errors, and the viewer should be updated accordingly. 
8. **Floating viewer** after applying layout:
   * Open new **demog** dataset.
   * Adjust Screen Zoom: Zoom out the screen.
   * Add and Move Viewer: Add a Trellis plot viewer. Detach the viewer and move it to the bottom of the screen.
   * Save the current layout.
   * Reset Screen Zoom: Zoom in the screen to its original size.
   * Apply the previously saved layout.
   * **_Expected Results:_**
     * The layout should automatically adjust to fit the current screen size.
     * All viewers should remain accessible and properly positioned on the screen.
9. **Split updates in properties**. 
   * Test on the SPGI dataset with opened new Trellis plot viewer: 
    * Open the Trellis plot's Context panel.
    * Add or remove X or Y column from the  Trellis plot itself.
    * Check the number of X or Y columns on the Context panel.
   * **_Expected results_**: number of selected columns in properties should be updated immediately on adding or removing a column.
10. **Pick Up / Apply**  ([#2116](https://github.com/datagrok-ai/public/issues/2116))
  * Open SPGI 
  * Add two trellis plots
  * Right-click the first trellis plot. Select **Pick up/Apply > Pick up** 
  * Right-click the second trellis plot. Select **Pick up/Apply > Apply**
  * Change the X or Y column on the first trellis plot. Observe whether these changes affect the second trellis plot
  * Adjust the zoom slider on the second trellis plot. Observe whether these changes affect the first trellis plot
  * **_Expected Result_**: After the initial settings transfer, the configurations of the two trellis plots should be completely independent. Changes to the X or Y column on the first plot or adjustments to the zoom slider on the second plot should not affect the other plot.
---
{
  "order": 13,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
