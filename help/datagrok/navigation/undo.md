---
title: "Undo and redo"
---

Press <kbd>Ctrl+Z</kbd> to reverse the last operation, and <kbd>Ctrl+Shift+Z</kbd> (or
<kbd>Ctrl+Y</kbd>) to re-apply it. The same commands are in the **Edit** menu. On macOS, use
<kbd>Cmd</kbd> instead of <kbd>Ctrl</kbd>.

Undo is **scoped to what you are looking at**: <kbd>Ctrl+Z</kbd> reverses the last operation on the
current table or view, never one you performed somewhere else. Switching to another view and back
brings the history with you.

## What can be undone

| Area | Operations |
|---|---|
| Grid editing | Editing a cell, clearing a cell (<kbd>Del</kbd>), toggling a checkbox, pasting a cell or a block |
| Rows | Removing selected rows, the grid's `−` and `+` icons, **Edit \| Add Rows** |
| Columns | Removing columns, changing a column's type, adding or editing a formula column, renaming a column and editing its properties |
| Values | **Find and Replace** |
| Metadata | Color coding and format changes from the column menu |
| Selection | **Reset Selection Filter** (<kbd>Esc</kbd>) |
| Layout | Closing a viewer or a view, adding a viewer, **View \| Layout \| Clear** |
| AI | Operations performed by the assistant (add/close viewer, change properties, filter, select) |

Everything else — sorting, filtering, creating new tables, and anonymizing data — is not recorded.
Operations that produce a new table are non-destructive, so nothing is lost by not undoing them.

## Editors with their own undo

Sketchers, the form designer, molecule and sequence editors, and graph editors keep their own undo
stacks. While the focus is inside one of them, <kbd>Ctrl+Z</kbd> goes to that editor; everywhere
else it goes to the platform. <kbd>Ctrl+Z</kbd> in a text field always does what it does in any
browser — it undoes typing, and it does not touch the table.

**While a dialog is open, <kbd>Ctrl+Z</kbd> never reaches the table.** A dialog is its own editing
session, and undoing something underneath it while you are filling in a form is never what you
meant. Close the dialog first.

## For developers

Plugins can make their own operations undoable:

```js
const before = column.get(0);
column.set(0, 42);
DG.UndoService.push('Set value', () => column.set(0, before), {
  redo: () => column.set(0, 42),
  context: column.dataFrame
});
```

- `grok.shell.undo()` / `grok.shell.redo()` — run the commands programmatically
- `grok.shell.canUndo` / `canRedo` — enablement
- `grok.shell.onUndo` / `onRedo` — observables carrying `{name, context}`
- `DG.UndoService.ownScope(element)` — claim <kbd>Ctrl+Z</kbd> for an embedded editor while focus is
  inside it, instead of installing a document-level key handler

Records hold values, not row indexes of live objects: capture what you need before the change, and
pass `context` so the record is released with its table.

See also: [Keyboard shortcuts](shortcuts.md).
