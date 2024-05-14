---
title: Controls
keywords:
 - shortcuts
 - controls
 - hot keys
sidebar position: 5
format: mdx
---

This article lists key Datagrok controls. Certain tools
or viewers may have additional shortcuts or modify the listed shortcuts for
different actions. For details, see the documentation for such tool or viewer.

* <kbd>Esc</kbd>:

    * **Viewers** | Clear selection _or_ toggle filters without resetting them. <br/> If
    you have both selected and filtered rows, pressing **Esc** first clears the
    selection, then toggles the filters.
     * **Full Screen mode** | Exit 


* <kbd>F1</kbd>: <b>Context Help</b> | Toggle
* <kbd>F2</kbd>: 

  * <b>Viewers | Grid</b> | Open column properties
  * <b>Viewers | Other than grid</b> | Edit title

* <kbd>F4</kbd>: <b>Context Panel</b> | Toggle
* <kbd>F5</kbd>: <b>Editors (e.g., scripts, queries)</b> | Run (requires an open editor)
* <kbd>F7</kbd>: <b>Presentation mode</b> | Toggle
* <kbd>F11</kbd>: <b>Full Screen mode</b> | Toggle  

<!-- Unassigned:
* <kbd>F6</kbd>: --
* <kbd>F3</kbd>: --
* <kbd>F8</kbd>: --
* <kbd>F9</kbd>: --
* <kbd>F10</kbd>: --
-->
<br/>

* <kbd>L</kbd>: <b>Viewers | Scatterplot</b> | Toggle rectangle/lasso selection mode
* <kbd>R</kbd>: <b>Viewers | Scatterplot</b> | Toggle regression line

<br/>

* <kbd>Ctrl + A</kbd>: <b>Viewers</b> | Select all
* <kbd>Ctrl + C</kbd>: <b>Viewers | Grid</b> | Copy cell
* <kbd>Ctrl + H</kbd>: <b>Viewers | Grid </b>| Open the <b>Find and replace</b> dialog 
* <kbd>Ctrl + L</kbd>: <b>Table View | Layouts</b> | Open Gallery (requires an open Table View)
* <kbd>Ctrl + O</kbd>: <b>Browse</b> | Open the <b>Open File</b> dialog
* <kbd>Ctrl + S</kbd>: <b>Table View | Layouts</b> | Save to Gallery (requires an open Table View)
* <kbd>Ctrl + V</kbd>: <b>Viewers | Grid</b> | Paste into cell
* <kbd>Ctrl + Z</kbd>: <b>Viewers | Grid</b> | Undo. You can reverse only one action. Supported actions:

   * Manually changing a cell value via a double-click
   * Clearing the [layout](../../visualize/view-layout.md)
   * Deleting rows or columns
   * Closing a viewer by clicking on the **x** icon
 
<!--
* <kbd>Ctrl + J</kbd>: Viewers | Grid | Undo. You can reverse only one action.
-->

* <kbd>Ctrl + Shift + A</kbd>: <b>Viewers</b> | Deselect all
* <kbd>Ctrl + Shift + R</kbd>: <b>Table View | Layouts</b> | Clear everything from a **Table View** except the grid
* <kbd>Ctrl + Shift + V</kbd>: <b>Table View</b> | Open a new view that is attached to the same table
* <kbd>Ctrl + Shift + Z</kbd>: <b>Viewers | Grid</b> | Redo. You can use the <b>Redo</b> (Ctrl+Shift+Z) command only after the <b>Undo</b> (Ctrl+Z) command

<br/>

* <kbd>Alt + C</kbd>: <b>Table View | Colum Manager</b> | Toggle
* <kbd>Alt + F</kbd>: <b>Full screen mode</b> | Open
* <kbd>Alt + S</kbd>: <b>Sidebar | Settings</b> | Open
* <kbd>Alt + T</kbd>: <b>Table Manager</b> | Toggle
* <kbd>Alt + V</kbd>: <b>Table View | Variables Panel</b> | Toggle
* <kbd>Alt + X</kbd>: <b>Table View | Toolbox</b> | Toggle

<!--
* Alt + Q: Search commands       Not working now
* Alt + I: Tools | Inspector | Open

TBD: Add Alt+ I +R/C to add new row or column 
-->

<br/>

* <kbd>Click</kbd>:

   * **Browse** | Make an object current (e.g., file or query)
   * **Viewers** | Make a row or cell object current 

* <kbd>Double click</kbd>: 

   * **Browse** | Open an object (e.g., file or query)
   * **Viewers | Grid** 
     * Cell | Edit a cell
     * Column header | Sort a column (Also works for grid-based viewers like a [heatmap](../../visualize/viewers/heat-map.md) or [correlation plot](../../visualize//viewers/correlation-plot.md))
   * **Viewers | Empty area** | Reset view in these viewers: 
   
     * Scatterplot
     * 3D Scatterplot
     * Histogram
     * Bar chart
     * Pie chart
     * Network diagram

* <kbd>Right-click</kbd>: Show context menu
* <kbd>Ctrl + Click</kbd>: <b>Viewers</b> | Toggle selected state, one at a time. For columns, Ctrl+click the column's header
* <kbd>Shift + Click</kbd>: <b>Viewers | Grid</b> | Select all rows above the one you click
* <kbd>Ctrl + Shift + Click</kbd>: <b>Viewers | Grid</b> | Clear selected rows above the row you click. For column headers, deselects columns

<br/>

* <kbd>Mouse drag</kbd>:

  * <b>Panel or View header</b> | Dock (detach, reposition)
  * <b>Viewers | Grid</b>:

    * Rows | Select rows and columns
    * Row headers | Select rows
    * Column header | Move a column. Selected columns
  are repositioned simultaneously next to each other
    * Right border of a row or column's header | Resize a row/column
  
  * <b>Viewers | Scatterplot, line chart</b>: Pan

* <kbd>Shift + Mouse drag</kbd>: <b>Viewers</b> | Select adjacent rows/data. For grid columns, Shift+drag column headers
* <kbd>Ctrl + Shift + Mouse drag</kbd>: <b>Viewers</b> | Deselect adjacent rows/data 
* <kbd>Alt + Mouse drag</kbd>: <b>Viewers | Scatterplot, line chart, bar chart</b> | Select an area to zoom in

<!--  
* <kbd>Ctrl + Shift + Mouse drag</kbd>: <b>Table Manager, Viewers (grid, scatterplot)</b> | ??? not working?
-->

<br/>

* <kbd>Mouse wheel up or down</kbd>:

  * <b>Views, Grid</b> | Scroll
  * <b>Viewers</b> | Zoom in/out for these viewers:
    
     * Scatterplot
     * 3D scatterplot
     * Line chart
     * Network diagram
     * Map


* <kbd>Up (↑), Down (↓)</kbd>: <b>Browse, Grid</b> | Navigate up/down
* <kbd>Left (←), Right (→)</kbd>:

  * <b>Browse</b>: Expand/collapse a tree node
  * <b>Viewers | Grid</b>: Navigate left/right

* <kbd>Ctrl + ↑↓</kbd>: <b>Viewers | Grid</b> | Go to previous/next selected row
* <kbd>Shift + ↑↓</kbd>: <b>Viewers | Grid</b> | Select consecutive rows
* <kbd>Ctrl + Shift + ↑↓</kbd>: <b>Viewers | Grid</b> | Deselect consecutive rows
* <kbd>PageUp, PageDown</kbd>: <b>Viewers | Grid</b> | Jump to the first/last row in the currently displayed row array. 
* <kbd>Home/End</kbd>: <b>Viewers | Grid</b> | Jump to the first/last column for the current row.
* <kbd>Ctrl + Home/End</kbd>: <b>Viewers | Grid</b> | Jump to the first/last row in the dataset.
* <kbd>Ctrl + Shift + Home/End</kbd>: <b>Viewers | Grid</b> | Select rows above/below current.

<br/>

* <kbd>Shift + Delete</kbd>: <b>Viewers | Grid</b> | Delete selected rows and/or columns
* <kbd>Ctrl + Enter</kbd>: <b>Editors (e.g., scripts, queries)</b> | Run (requires an open editor)

<br/>

:::tip tips

Some tools, like the [Table Manager](panels/table-manager.md) or 
[Column Manager](panels/column-manager.md), are based on the
[grid](../../visualize/viewers/grid.md). Consequently, grid shortcuts 
also apply to the corresponding elements within these tools (e.g., rows).

:::