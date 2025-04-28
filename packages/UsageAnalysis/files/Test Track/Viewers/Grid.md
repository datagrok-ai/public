### Grid

#### 1. Viewer basics 

1. Open SPGI
1. Double-click any column's header to sort it (check molecular, numeric, bool and string columns)
1. Resize some columns
1. Resize rows - **columns width should be adjusted automatically**
1. Select a row, then hold Shift or Ctrl to select multiple rows
1. Double-click a cell, change its value, and press Enter (check molecular, numeric, and string columns) - **the value is updated in the grid immediately**


#### 2. Columns

1. Open the hamburger menu
1. **Filter**: 
  * Use the filter 
  * Click Add filter - **the filter should be added to the Filter Panel**
1. **Color**: 
  * Turn on the color coding
  * Change the color coding type (where aplicable) 
  * Apply coloring to the text/background
1. Check the other menu items.
1. Repeate for molecular, numeric, bool and string columns

#### 4. Column's context menu 

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

#### 3. Context Panel

 **Column**
1. Open the **Context Panel** 
1. Click any column's header
1. Open the Filter Panel
1. Check the sync between filters on the Context and Filter Panels
1. Check color coding sync between the Context Panel and the column's menus
1. Check the ability to change style
1. Go to **Advanced** > **Permissions** > **Edited by** and enter your user name - **you should have the ability to edit the column**
1. Go to Advanced > Permissions > Edited by and enter only your colleague's user name - **you should NOT have the ability to edit the column, a balloon should appear**
 1. Check the other menu items

**Grid**
1. Click the Gear icon
1. Modify the settings and check that the changes are reflected in the grid

#### 4. Context menu

1. Right-click the grid
1. Use **Add** item and add 
  * Column
  * Summary column (all of them)
  * Top
  * Column Stats
1. Go to the **Chemical Space X** and rename it - check previously added columns
1. Use **Context menu > Add > Summary column > Smartform**
1. Edit the smartform and check its properties on the Context Panel
1. Open SPGI-linked1. Use Data > Link tables and link SPGI and SPGI-linked1.
1. Go to SPGI
1. Use **Context menu > Add > Linked Tables > SPGI-linked1** - check the result
1. Open the **Order Or Hide Columns..** dialog and check its functionality

#### 6. Pick Up / Apply   ([#1887](https://github.com/datagrok-ai/public/issues/1887))

1. For **Average mass** set **Color Coding** and edit color scheme
3. **Pick Up/ Apply > Pick Up**
3. For **TPSA** column, use **Pick Up/ Apply > Apply** - The **TPSA** column should match the **Average mass** column’s color-coding 
4. Right-click the grid and select **Grid Color Coding > All**
1. Go to the **Context Panel > Misc > Color Scheme** and change it
1. **Pick Up/ Apply > Pick Up** 
5. Open the second SPGI dataset
1. Right-click the grid and select **Pick Up/Apply > Apply**. The color-coding, formatting and style on both grids should match.

#### 7. Column groups

1. Select some columns 
1. Remove them
1. Change their properties

#### 8. Filtering

1. Open the **Filter Panel**  
2. Apply structure, categorical, and numeric filters.  
**Verify: The grid updates continuously as filters are adjusted**

#### 7. Layout and Project save/restore
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