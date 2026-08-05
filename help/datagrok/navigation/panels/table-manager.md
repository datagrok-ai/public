---
title: "Table Manager"
unlisted: true
sidebar_position: 3
description: Use Table Manager to view, select, and batch-manage all currently open tables in a grid.
keywords:
  - view tables
  - manage open tables
  - table list
  - batch actions on tables
  - multiple table selection
  - table metadata
---

Table Manager, available via **View | Tables**, contains a list of currently open tables. Use it to navigate between
tables, select them, or perform batch actions. You can also view metadata for multiple tables in a tabular format.

The implementation is based on the [grid](../../../visualize/viewers/grid.md), so many of the grid's features apply.

Usage:

|                      |                |
|----------------------|----------------|
| Click                | Activate or create table view   |
| Shift+drag           | Select multiple tables |
| Ctrl+click           | Toggle table selection |
| Esc                  | Clear table selection |
| Right-click          | Show popup menu |
| Popup: Show \        | Toggle visibility of the property "x" |
| Popup: Open as table | Add tables to workspace (tables in rows) |

If a context menu is open when multiple tables are selected, you can choose to apply commands to
either the current table or all selected tables.

See also:

* [Column Manager](column-manager.md)
* [Grid](../../../visualize/viewers/grid.md)
