<!-- TITLE: Tests: Viewers -->
<!-- SUBTITLE: -->

# Tests: Viewers

A [viewer](../viewers.md) is a visual component associated with a table.

List of viewers available for testing:

* [Scatter Plot](../visualize/viewers/scatter-plot.md)
* [Histogram](../visualize/viewers/shistogram.md)
* [Line Chart](../visualize/viewers/sline-chart.md)
* [Bar chart](../visualize/viewers/sbar-chart.md)
* [Box plot](../visualize/viewers/sbox-plot.md)
* [Filters](../visualize/viewers/sfilters.md)
* [Trellis Plot](../visualize/viewers/strellis-plot.md)
* [Tree Map](../visualize/viewers/stree-map.md)
* [Calendar](../visualize/viewers/scalendar.md)
* [Google Map](../visualize/viewers/sgoogle-map.md)
* [Grid](../visualize/viewers/sgrid.md)
* [3D Scatter Plot](../visualize/viewers/3d-scatter-plot.md)
* [Matrix Plot](../visualize/viewers/smatrix-plot.md)
* [Network Diagram](../visualize/viewers/snetwork-diagram.md)
* [Parallel Coordinates Plot](../visualize/viewers/spc-plot.md)
* [Pie Chart](../visualize/viewers/spie-chart.md)
* [Scripting Viewer](../visualize/viewers/scripting-viewer.md)
* [Statistics](../visualize/viewers/sstatistics.md)
* [Markup Viewer](../visualize/viewers/smarkup-viewer.md)
* [Word Cloud](../visualize/viewers/sword-cloud.md)
* [Scripting Viewer](../visualize/viewers/scripting-viewer.md)

Approaches to testing, Testing scenarios, test items in general are similar for all viewers. Categories of testing scenarios common to all viewers were identified:

1. Creating a [viewer](../viewers.md). (UX, tooltips, graphic icons, on-screen display):

   * From toolbox
   * From menu "Add"

1. Correct construction of created [viewer](../viewers.md). 
Here we look at the correspondence of the chosen viewer type, correctness of construction from the point of view of math logic.

1. Working with a [viewer](../viewers.md) window:

   * Moving the window by the workspace
   * Zoom in\Zoom out
   * Navigation inside viewer window
   * Hover and select elements inside viewer
   * Open\Close property window
   * Display of help window
   * The work of the controls located in the viewer area (e.g. popup menu "Size" on Scatter Plot)
   * Return viewer after closing by **Edit | Undo** (```CTRL + Z```)
   
1. Property window
   * Appearance and UI
   * Moving the property window by the workspace
   * Property search
   * Change the property values of the "Data" item
   * Change the property values of the "Appearance" item
   * For [Scripting Viewer](../visualize/viewers/scripting-viewer.md):
       * Edit script in **Appearance | Script** 

1. Context menu
   * Call the context menu
   * Navigation through the items of context menu
   * Change viewer settings from the context menu

1. "Viewer" submenu
   
   * Copy markup
   * Properties
   * Clone
   * Use in Trellis
   * Save to Gallery
   * Embed
   * Full Screen
   * Close

1. "Style" submenu
   * Pick up style
   * Apply data settings
   * Apply style settings
   * Apply
   * Set as default
   * Reset default

1. "Tooltip" submenu
   * Use a group tooltip
   * Reset 

1. Used input data

   * Different types  (string, int, double, datetime, percent, bool)
   * Full and non-full (with nulls or empty tables) data
   * Selected rows should be highlighted on the viewer

## Sample test-case for testing viewers on example of "Scatter Plot"

1. Open *"demog"* table

1. Add "Scatter plot" viewer to layout (from "Viewers" tab on Toolbax or from "Add" menu)
   * "Scatter Plot" added on layout
   * Help switched to "Scatter Plot" page
   * *"Weight"* column for X axis and *"Height"* column for Y axis are automatically established
   
1. Select *"Age"* column in "Size" selector on [Scatter Plot](../visualize/viewers/scatter-plot.md)
   * Points on [Scatter Plot](../visualize/viewers/scatter-plot.md) displayed in different sizes depending on values ​​of *"Age"* column
   * Under "Size" selector only columns with numerical types are represented
   
1. Drag *"Sex"* column from grid to "Color" selector on [Scatter Plot](../visualize/viewers/scatter-plot.md)
   * Points on [Scatter Plot](../visualize/viewers/scatter-plot.md) displayed in different sizes depending on values ​​of *"Age"* column
   * When dragging column, "Colort" selector on [Scatter Plot](../visualize/viewers/scatter-plot.md) is highlighted

1. Delete *"Height" column" form *"demog"* table
   * Y axis on [Scatter Plot](../visualize/viewers/scatter-plot.md) switched to another suitable column

1. Delete *"Sex"* column" form *"demog"* table 
   * Color mapping off on [Scatter Plot](../visualize/viewers/scatter-plot.md)
   
1. Open properties for [Scatter Plot](../visualize/viewers/scatter-plot.md)
   *  Viewer properties are displayed on [Property Panel](../../overview/navigation.md#properties)
   
1. Check performance and impact of changing properties on [Scatter Plot](../visualize/viewers/scatter-plot.md)
   * Data and appearance display properties changing the for  [Scatter Plot](../visualize/viewers/scatter-plot.md) affects expected
   * After changing appearance properties the [Scatter Plot](../visualize/viewers/scatter-plot.md) does not displayed "default"
   
1. Pick up customized style of [Scatter Plot](../visualize/viewers/scatter-plot.md) as default from **Context menu | Style | Set as default**

1. Add one more [Scatter Plot](../visualize/viewers/scatter-plot.md) on layout
   * New [Scatter Plot](../visualize/viewers/scatter-plot.md) has data and appearance settings as configured in step 8

1. Delete [Scatter Plot](../visualize/viewers/scatter-plot.md) from previous step

1. Reset view of [Scatter Plot](../visualize/viewers/scatter-plot.md) from it's context menu
   * [Scatter Plot](../visualize/viewers/scatter-plot.md) style is now "default"
   
1. Use  mouse wheel to zoom [Scatter Plot](../visualize/viewers/scatter-plot.md) area
   * When scrolling wheel "away from you", [Scatter Plot](../visualize/viewers/scatter-plot.md) area zooming to current padding of the cursor
   
1. Select [Scatter Plot](../visualize/viewers/scatter-plot.md) area with ```Alt``` pressed   
   * Selected area is displayed on [Scatter Plot](../visualize/viewers/scatter-plot.md) 
   
1. Double click on [Scatter Plot](../visualize/viewers/scatter-plot.md) area   
   * [Scatter Plot](../visualize/viewers/scatter-plot.md) area  scale returned to default
   
1. Use [Scatter Plot](../visualize/viewers/scatter-plot.md) for tooltip of other viewers. Select **Context menu | Tooltip | Use as group tooltip**  

1. Add [Histogram](../visualize/viewers/histogram.md) viewer to layout
   * When hovering over histogram group tooltip is displayed, which shows [Scatter Plot](../visualize/viewers/scatter-plot.md) that corresponds to group
   
1. Delete all columns from *"demog"* table 
   * Nothing is shown on [Scatter Plot](../visualize/viewers/scatter-plot.md) view  
   * No errors and crashes 

See also:
 * [Viewers](../viewers.md)
 * [Creating dashboards tutorial](../tutorials/creating-dashboards.md)
