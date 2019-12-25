<!-- TITLE: Table Manager -->
<!-- SUBTITLE: -->

# Table Manager

Table Manager, available via **View | Tables**, contains a list of currently open tables. Use it to 
navigate between tables, select them, or perform batch actions. It also allows to view metadata on
multiple tables in a tabular format.

The implementation is based on the [grid](../viewers/grid.md), so many of the grid's features apply.

Usage:

|                  |                |
|------------------|----------------|
| Click            | Activate or create table view   |
| Shift+drag       | Select multiple tables |
| Ctrl+click       | Toggle table selection |
| Esc              | Clear table selection |
| Right-click      | Show popup menu |
| Popup: Show \| x | Toggle visibility of the property "x" |
| Popup: Save as table | Add tables to workspace (tables in rows) |

If a context menu is open when multiple columns are selected, user will be present with a 
choice to apply commands to either current table, or all selected tables.

See also:
* [Column Manager](column-manager.md)
* [Grid](../viewers/grid.md)
