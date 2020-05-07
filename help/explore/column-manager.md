<!-- TITLE: Column Manager -->
<!-- SUBTITLE: -->

# Column Manager

Column Manager, available via **View | Columns** or **Add | Column Viewer**, contains a list of 
columns in the currently open tables. Use it to navigate between columns, select them, or perform batch actions. 
It also allows to view metadata on multiple columns in a tabular format, as well as other column
properties, such as statistics.

The implementation is based on the [grid](../visualize/viewers/grid.md), so many of the grid's features apply.

Usage:

|                  |                |
|------------------|----------------|
| Click            | Jump to column |
| Shift+drag       | Select multiple columns |
| Ctrl+click       | Toggle column selection |
| Esc              | Clear selection |
| Right-click      | Show popup menu |
| Popup: Add Stats | Show/hide statistics |

If a context menu is open when multiple columns are selected, user will be present with a 
choice to apply commands to either current table, or all selected tables.

See also:

* [Grid](../visualize/viewers/grid.md)
