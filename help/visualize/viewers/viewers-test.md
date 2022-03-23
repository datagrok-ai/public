<!-- TITLE: Tests: Viewers -->
<!-- SUBTITLE: -->

# Tests: Viewers

A [viewer](../viewers.md) is a visual component associated with a table.

List of viewers available for testing:

* [Scatter plot](../viewers/scatter-plot.md)
* [Histogram](../viewers/histogram.md)
* [Line chart](../viewers/line-chart.md)
* [Bar chart](../viewers/bar-chart.md)
* [Box plot](../viewers/box-plot.md)
* [Filters](../viewers/filters.md)
* [Trellis plot](../viewers/trellis-plot.md)
* [Tree map](../viewers/tree-map.md)
* [Calendar](../viewers/calendar.md)
* [Google map](../viewers/google-map.md)
* [Grid](../viewers/grid.md)
* [3D Scatter plot](../viewers/3d-scatter-plot.md)
* [Matrix plot](../viewers/matrix-plot.md)
* [Network diagram](../viewers/network-diagram.md)
* [Parallel coordinates plot](../viewers/pc-plot.md)
* [Pie chart](../viewers/pie-chart.md)
* [Scripting viewer](../viewers/scripting-viewer.md)
* [Statistics](../viewers/statistics.md)
* [Word cloud](../viewers/word-cloud.md)
* [Scripting viewer](../viewers/scripting-viewer.md)

Approaches to testing, Testing scenarios, test items in general are similar for all viewers. Categories of testing
scenarios common to all viewers were identified:

1. Creating a [viewer](../viewers.md). (UX, tooltips, graphic icons, on-screen display):

    * From toolbox
    * From menu "Add"

2. Correct construction of created [viewer](../viewers.md). Here we look at the correspondence of the chosen viewer
   type, correctness of construction from the point of view of math logic.

3. Working with a [viewer](../viewers.md) window:

    * Moving the window by the workspace
    * Zoom in\Zoom out
    * Navigation inside viewer window
    * Hover and select elements inside viewer
    * Open\Close property window
    * Display of help window
    * The work of the controls located in the viewer area (e.g. popup menu "Size" on Scatter Plot)
    * Return viewer after closing by **Edit | Undo** (```CTRL + Z```)

4. Property window

    * Appearance and UI
    * Moving the property window by the workspace
    * Property search
    * Change the property values of the "Data" item
    * Change the property values of the "Appearance" item
    * For [Scripting Viewer](../viewers/scripting-viewer.md):
    * Edit script in **Appearance | Script**

5. Context menu

    * Call the context menu
    * Navigation through the items of context menu
    * Change viewer settings from the context menu

6. "Viewer" submenu

    * Copy markup
    * Properties
    * Clone
    * Use in Trellis
    * Save to Gallery
    * Embed
    * Full Screen
    * Close

7. "Style" submenu

    * Pick up style
    * Apply data settings
    * Apply style settings
    * Apply
    * Set as default
    * Reset default

8. "Tooltip" submenu

    * Use a group tooltip
    * Reset

9. Used input data

    * Different types  (string, int, double, datetime, percent, bool)
    * Full and non-full (with nulls or empty tables) data
    * Selected rows should be highlighted on the viewer

## Sample test-case for testing viewers on example of "scatter plot"

1. Open *"demog"* table

2. Add "Scatter plot" viewer to layout (from "Viewers" tab on Toolbax or from "Add" menu)

    * "Scatter Plot" added on layout
    * Help switched to "Scatter Plot" page
    * *"Weight"* column for X axis and *"Height"* column for Y axis are automatically established

3. Select *"Age"* column in "Size" selector on [Scatter Plot](../viewers/scatter-plot.md)

    * Points on [Scatter Plot](../viewers/scatter-plot.md) displayed in different sizes depending on values ​​of *"Age"*
      column
    * Under "Size" selector only columns with numerical types are represented

4. Drag *"Sex"* column from grid to "Color" selector on [Scatter Plot](../viewers/scatter-plot.md)

    * Points on [Scatter Plot](../viewers/scatter-plot.md) displayed in different sizes depending on values ​​of *"Age"*
      column
    * When dragging column, "Colort" selector on [Scatter Plot](../viewers/scatter-plot.md)
      is highlighted

5. Delete *"Height"* column" form *"demog"* table

    * Y axis on [Scatter Plot](../viewers/scatter-plot.md) switched to another suitable column

6. Delete *"Sex"* column" form *"demog"* table

    * Color mapping off on [Scatter Plot](../viewers/scatter-plot.md)

7. Open properties for [Scatter Plot](../viewers/scatter-plot.md)

    * Viewer properties are displayed on [Property Panel](../../overview/navigation.md#properties)

8. Check performance and impact of changing properties on [Scatter Plot](../viewers/scatter-plot.md)

    * Data and appearance display properties changing the for  [Scatter Plot](../viewers/scatter-plot.md) affects
      expected
    * After changing appearance properties the [Scatter Plot](../viewers/scatter-plot.md)
      does not displayed "default"

9. Pick up customized style of [Scatter Plot](../viewers/scatter-plot.md) as default from **Context menu | Style | Set
   as default**

10. Add one more [Scatter Plot](../viewers/scatter-plot.md) on layout

    * New [Scatter Plot](../viewers/scatter-plot.md) has data and appearance settings as configured in step 8

11. Delete [Scatter Plot](../viewers/scatter-plot.md) from previous step

12. Reset view of [Scatter Plot](../viewers/scatter-plot.md) from it's context menu

    * [Scatter Plot](../viewers/scatter-plot.md) style is now "default"

13. Use mouse wheel to zoom [Scatter Plot](../viewers/scatter-plot.md) area

    * When scrolling wheel "away from you", [Scatter Plot](../viewers/scatter-plot.md) area zooming to current padding
      of the cursor

14. Select [Scatter Plot](../viewers/scatter-plot.md) area with ```Alt``` pressed

    * Selected area is displayed on [Scatter Plot](../viewers/scatter-plot.md)

15. Double click on [Scatter Plot](../viewers/scatter-plot.md) area

    * [Scatter Plot](../viewers/scatter-plot.md) area scale returned to default

16. Use [Scatter Plot](../viewers/scatter-plot.md) for tooltip of other viewers. Select **
    Context menu | Tooltip | Use as group tooltip**

17. Add [Histogram](../viewers/histogram.md) viewer to layout

    * When hovering over histogram group tooltip is displayed, which shows [Scatter Plot](../viewers/scatter-plot.md)
      that corresponds to group

18. Delete all columns from *"demog"* table

    * Nothing is shown on [Scatter Plot](../viewers/scatter-plot.md) view
    * No errors and crashes

See also:

* [Viewers](../viewers.md)
