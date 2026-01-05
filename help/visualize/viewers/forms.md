---
title: "Forms"
---

**Forms viewer** shows rows in forms positioned side-by-side.
It is useful for quick exploration of the selection, as well as for comparison between current, 
mouse-over, and selected rows.

:::note Note

To use a forms viewer, install the [PowerGrid](https://github.com/datagrok-ai/public/tree/master/packages/PowerGrid) package.

:::

* To add **Forms**, click the "Add viewer" icon on top, and then select "Forms".
* To select columns on a form, use the **Fields** property in the options. Drag-and-drop columns to reorder. 
* To choose a subset of rows to show, use these properties: `Show Current Row`, `Show Mouse Over Row`, and 
  `Show Selected Rows`. You can mix and match these options. 

Green color stripe on top of the form indicates current row, grey one indicates mouse-over row. 
The viewer also works together with the grid:

* To select or deselect rows, Ctrl+click the form.
* To go to a particular cell in a grid, click on a field
* To make a column current, click on the column name.

![Forms viewer](img/forms.gif)

To reorder the fields in Forms viewer, go to Froms viewer settings, open Fileds option and just drag-n-drop corresponding fields in the ‘Select columns…’ dialog.

![Forms viewer](img/reorder_fields_in_forms_viewer.gif)

## Sorting

The forms viewer displays rows in the following priority order 
(if the corresponding options are enabled):

1. Current row  
2. Mouse-over row  
3. Selected rows

By default, the order of selected rows is inherited from the associated grid. 
Any changes to the row order in the grid are automatically reflected in the forms viewer.

To sort the selected rows by a column directly in the forms viewer, double-click the column header 
or use the **Sort by** setting (**Context Panel > Data**). A “↓” symbol appears next to the column name 
to indicate sorting in descending order. Double-click again to sort ascending. The next double-click resets sorting.

![Forms viewer sorting](img/forms-sorting.gif)

See also:
* [Viewers](../viewers/viewers.md)
* [Table View](../table-view-1.md)
* [Community: Visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521)
