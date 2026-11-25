---
name: datagrok-undo
description: How to change an open table so that Ctrl+Z (Edit | Undo) takes the change back, and what undo cannot do. Use before adding, removing or renaming columns, removing rows, overwriting a column's values, changing column metadata (format, color coding, tags), adding or closing viewers, clearing the layout or closing a view, and when the user asks to undo or redo something or asks what can be undone.
---

# datagrok-undo

The platform undoes its own menu actions. A change a script makes directly (`t.columns.remove`,
`t.rows.removeWhereIdx`, `col.init`, `view.addViewer`) is **not** undoable. Three ways make one that is, in this
order of preference.

## 1. The table tools

Each acts on the table of the current table view, records one undo step and returns what it read back, with
`undoable: true | false`.

| Change | Tool |
|---|---|
| remove columns | `table_remove_columns {columns}` |
| remove rows | `table_remove_rows {where}` |
| add a column from a per-row value or a formula | `table_add_column {name, type, value}` / `{name, formula}` |
| overwrite a column's values | `table_set_values {column, value}` |
| rename; format, units, semantic type, color coding, tags | `table_edit_column {column, ...}` |
| add a viewer | `table_add_viewer {type, options}` |

`where` and `value` are JS expressions run for every row: `i` is the row index, `t` the table, `v` the current
value (`null` when missing). Like any action, a tool call is followed by `datagrok_verify`.

## 2. The `undoable` global of `datagrok_exec`

For a change a script has to make itself, e.g. a column a platform function added:

| Change | Call |
|---|---|
| a column the script or a function it called has added | `undoable.columnAdded(col)` |
| anything that rewrites a column's values | `const ok = await undoable.changeValues(col, async () => { ... })` |
| anything that renames a column or sets its tags | `await undoable.editColumn(col, () => { ... })` |
| remove columns / rows, add a viewer | `undoable.removeColumns(t, names)`, `await undoable.removeRows(t, (i) => ...)`, `undoable.addViewer(view, type, options)` |

## 3. Platform actions that record their own step

| Change | Call |
|---|---|
| clear the row and column selection | `await grok.functions.call('CmdResetSelectionFilter')` |
| remove every viewer but the grid | `grok.shell.topMenu.find('View').find('Layout').find('Clear').click()` |
| close a table view | `view.close()` |

## Limits of undo

- **Not undoable at all:** filtering, sorting, a viewer's properties, closing one viewer (`v.close()`),
  `view.resetLayout()`, opening a table, joins and aggregations (they make new tables), anything saved on the
  server (projects, files, queries, shares), and every change made around the three ways above, including one
  made through `call_view_function`.
- **One step per tool or helper call.** A request that made three changes takes three Ctrl+Z, the latest first.
  A step in the middle cannot be undone on its own.
- **Size.** Removed rows are kept only up to 20 MB: above that the platform says "too large to undo", removes
  them for good and forgets the table's earlier steps. Columns are kept or copied only up to 2.5 million cells
  (rows × columns): beyond that they are removed or overwritten without a step. Each case returns
  `undoable: false`.
- **Nothing changed.** A call that changed nothing (`removed: 0`, an edit to the same values) records no step
  and also returns `undoable: false`.
- **History.** The last 20 steps are kept, 50 MB in all, in the browser's memory only: reloading the page or
  closing the table loses them. A new change clears redo. Overwritten values have undo, not redo.
- **Stale steps.** Steps belong to a table. If an unrecorded change altered its rows since, the platform skips
  the older steps ("data has changed").
- **Selection.** Undoing a row removal brings the rows back selected.
- **Keyboard.** Ctrl+Z (Cmd+Z on a Mac) in the AI prompt box edits the prompt text. The table needs focus: a
  click on it first, or `Edit | Undo`.

## What to tell the user

Nothing: do not mention undo in a reply, neither that a change can be undone nor that it cannot. The exception
is a user who asks about it:

- "Can I undo this?" / "how do I take it back?": answer from the limits above; the shortcut is "Ctrl+Z (Cmd+Z on
  a Mac), after a click on the table, or Edit | Undo".
- "Undo that" / "redo": read `DG.UndoService.undoName` (`redoName`). If it names the change the user means, call
  `grok.shell.undo()` (`redo()`) and say what was undone; otherwise say that this change cannot be undone.
