---
title: Ten thousand rows
group: Lists & tables
order: 30
summary: A virtualized table with a sticky header, sortable columns and keyboard selection.
live: yes
---

`DataTable` is for more rows than the DOM should hold: only the visible window is rendered, cells pooled with
their row. Columns are data — a `name` (the item's property), a header, a CSS grid width, a `render` — and a
header click writes the `sort` signal (`'ic50'`, `'!ic50'` descending, `''` none); the rows are a `computed` over
it, so the order is yours to decide. Double-click or Enter activates a row; a hovered row is tinted. Give flexible
columns a floor (`width: 'minmax(90px, 1fr)'`): in a narrow pane the table then scrolls sideways instead of
squeezing a column to nothing. A cell without `render` shows `String(value)`, so the Measured column draws its date with
`Format.date`: the user's platform date format, read in the browser's zone.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (the height, on the table's own root). u2 handles the
virtualization, the sticky header, the sort arrows and `aria-sort`, column resizing, selection with the arrows,
Page Up/Down and Home/End or with Ctrl- and Shift-click, and activation.
