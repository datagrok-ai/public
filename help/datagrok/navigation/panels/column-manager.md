---
title: "Column Manager"
unlisted: true
sidebar_position: 4
---

Column Manager, available via **View | Columns** or **Add | Column Viewer**, contains a list of columns in the currently
open tables. Use it to navigate between columns, select them, or perform batch actions. You can also view metadata
for multiple columns in a tabular format, as well as other column properties, such as statistics.

The implementation is based on the [grid](../../../visualize/viewers/grid.md), so many of the grid's features apply.

Usage:

|                  |                |
|------------------|----------------|
| Click            | Jump to column |
| Shift+drag       | Select multiple columns |
| Ctrl+click       | Toggle column selection |
| Esc              | Clear selection |
| Right-click      | Show popup menu |
| Popup: Add Stats | Show/hide statistics |

If a context menu is open when multiple columns are selected, you can choose to apply commands to
either the current table or all selected tables.

See also:

* [Grid](../../../visualize/viewers/grid.md)
