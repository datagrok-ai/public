---
title: A panel bound to the table's current row
group: Inside Datagrok
order: 20
summary: A docked editor whose fields are the current row's cells — move the row, the fields follow; edit, the cell changes.
live: no
---

A spec binds straight to a workspace table: a `u2-table-source` node names it, and `$.table.currentRow.AGE` is
that cell of the current row — two-way, so typing writes the cell, and moving the current row in the grid (or
typing a row number into the first field) loads the next one. `$.table.selection`, `filter`, `rowCount` and
`columns` are there too. The editor is built under the view's host and docked with `host.dock(editor, 'right')`:
closing the panel, or the view, disposes the editor and every subscription it made, and disposing the editor
closes the panel.

## What it takes

{lines} lines for a live row editor docked beside the grid. u2 handles the current-row and cell signals (they
follow edits made anywhere), writing back, a closed table (the fields go empty instead of breaking), and the
disposal with the panel and the view.
