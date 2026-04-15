### Grid

#### 1. Viewer basics 

1. Open SPGI
1. Double-click any column's header to sort it (check molecular, numeric, bool and string columns) - **sort cycles asc/desc/unsort**
1. Resize some columns by dragging the border between column headers
1. Resize rows by dragging the bottom border of a row number - **row height changes for all rows**
1. Click a row, then Shift+Click another row to select a range; Ctrl+Click to toggle individual rows
1. Double-click a cell, change its value, and press Enter (check molecular, numeric, and string columns) - **the value is updated in the grid immediately**


#### 2. Columns

1. Hover a column header and click the hamburger (three-lines) icon to open the column menu
1. **Filter**:
  * Adjust the inline filter - **rows should be filtered immediately**
  * Click **Add filter** - **the filter is added to the Filter Panel**
1. **Color**:
  * Turn on color coding
  * Switch the color coding type (Linear / Categorical / Conditional, where applicable)
  * Toggle coloring between text and background
1. Check the other menu items
1. Repeat for molecular, numeric, bool and string columns

#### 3. Column's context menu 

1. Right-click the **Structure** column's header
1. Switch **Renderer** to `As Text` and back to `As Structure`
1. Right-click the **Series** column's header and check:
  * Sorting
  * Color coding (turn on/off, edit)
  * Hide (return visibility via the **Order Or Hide Columns..** dialog)
1. Right-click the **Chemical Space X** column's header and check:
  * Set a different format
  * Sort
  * Change type

#### 4. Context Panel

 **Column**
1. Open the **Context Panel** 
1. Click any column's header
1. Open the **Filter Panel**
1. Check the sync between filters on the **Context** and **Filter Panels**
1. Check color coding sync between the Context Panel and the column's menus
1. Check the ability to change style
1. Go to **Advanced** > **Permissions** > **Edited by** and enter your user name - **you should have the ability to edit the column**
1. Go to **Advanced** > **Permissions** > **Edited by** and enter only your colleague's user name - **you should NOT have the ability to edit the column, a balloon should appear when trying to edit**
 1. Check the other menu items

**Grid**
1. Click the **Gear** icon
1. Modify the settings and check that the changes are reflected in the grid

#### 5. Context menu

1. Right-click the grid
1. Use **Add** submenu and add:
  * **Column...**
  * **Summary Columns** (try each renderer)
  * **Top > Histogram** (and other Top options)
1. Rename the **Chemical Space X** column - **previously added summary columns still reference it correctly**
1. Use **Context menu > Add > Summary Columns > Smart Form**
1. Edit the Smart Form and check its properties on the Context Panel
1. Open SPGI-linked1, then use **Data > Link Tables...** to link SPGI and SPGI-linked1
1. Go back to SPGI
1. Use **Context menu > Add > Linked Tables > SPGI-linked1** - **linked columns appear in the grid**
1. Open **Order or Hide Columns...** dialog - toggle visibility, reorder, and apply

#### 6. Pick Up / Apply   ([#1887](https://github.com/datagrok-ai/public/issues/1887))

1. For **Average mass** set the linear **Color Coding** and edit color scheme
3. **Color Coding** >**Pick Up Coloring**
3. For **TPSA** column, use ***Color Coding** >**Apply Coloring** - The **TPSA** column should match the **Average mass** column’s color-coding 
4. Right-click the grid and select **Grid Color Coding > All**
1. Go to the **Context Panel > Misc > Color Scheme** and change it
1. **Pick Up/ Apply > Pick Up** 
5. Open the second SPGI dataset
1. Right-click the grid and select **Pick Up/Apply > Apply**. The color-coding, formatting and style on both grids should match.

#### 7. Column groups

1. Right-click a column header and choose **Group Columns...** (or select several columns and group them)
1. Expand/collapse the group in the header
1. Remove the group via right-click > **Ungroup**
1. Change group properties (name, color) via the Context Panel

#### 8. Filtering

1. Open the **Filter Panel**  
2. Apply structure, categorical, and numeric filters

**Verify: The grid updates continuously as filters are adjusted**

#### 9. Layout and Project save/restore
1. Save the layout
2. Change the layout (e.g. add some viewers)
1. Apply the saved layout - **Verify that all formatting, styling, coloring, column order, and other changes are preserved.**
1. Save the project
1. Close All
1. Open the saved project - **check**

---
{
  "order": 25,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}