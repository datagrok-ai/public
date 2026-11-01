---
title: Rows the user edits
group: Inputs & validation
order: 40
summary: A template's fields as editable, reorderable rows — each cell a real input, the value the records.
live: yes
---

`RecordsInput` is the input for a short list of typed rows: a column per field, each cell an inline input from
the column's `input`, add, remove and — with `orderable` — drag or keyboard reordering. Its value is the records,
one new array per edit with the untouched records kept as they were, so it validates and saves like a text box.
A cell's verdict is its red underline; the input's own line names the row ("Field 3: Value can't be empty").

## What it takes

{lines} lines. u2 handles the row editors, add/remove/reorder (with Alt+arrows on the keyboard), focus moving to a
new row, per-cell and whole-list validation, and the reason on Save.
