---
title: A dialog with a form
group: Dialogs
order: 10
summary: Required fields gate OK with the reason; an async OK keeps the dialog busy and a failure stays in it.
live: yes
height: 380
old: ui.dialog today
---

A `Dialog` with a `Form` in it waits for the form: OK is disabled while a field fails, and its tooltip names the
field ("Name: Value can't be empty"). An `onOK` that returns a promise keeps the dialog open and busy — the footer
disabled, a spinner on OK — and closes it when the promise resolves; a rejection stays in the dialog's error row,
with the user's input intact. Escape or ✕ while it runs aborts `signal`. CANCEL comes with the first action, as in
every platform dialog; a dialog whose work applies as it is done takes `{closeOnly: true}` and has one CLOSE
instead. A dialog is a root: built by a button's handler, it outlives the button's row. Register "Aspirin" to see
the failure.

## What it takes

{lines} lines against {oldLines} with `ui.dialog` (second tab), which wires the busy state, the error and the
button's off state by hand. u2 handles the focus (the first field on open, back to the button on close),
Enter for OK and Escape for CANCEL, dragging by the title, and disposing the inputs on close.
