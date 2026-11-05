---
title: A pick with an order
group: Inputs & validation
order: 50
summary: Which steps to run and in what order — one checklist whose value is the checked items, in order.
live: yes
---

`ChecklistInput({orderable: true})` is for a subset where the order matters too: the functions to run, the
columns to show. Its `value` is the checked items in the order shown; `order` is every item in that order, so
what is left out is `items` minus `value`. Drag a row by its handle or press Alt+↑/↓ in it.

## What it takes

{lines} lines. u2 handles the checkboxes, the drag handle, keyboard reordering with focus kept on the moved row,
and `required` (at least one checked).
