---
title: A list with details
group: Panels & layout
order: 30
summary: Forty compounds beside the current one's details, split by a sash: a MasterDetail.
live: yes
old: ui.splitH today
---

A list beside its current item's details is a `MasterDetail`: `row` says what each row shows (a label, and a `hint`
after it), `detail` builds the current item's pane the first time it is shown and keeps it, and the user drags the
sash or moves the current row with ↑/↓. Its height is `--u2-master-detail-height`, set on an element of the
sample's own; it fills its host by default. A detail pane on its own is a `tableFromMap`: muted keys, values that
wrap, and no hover tint — it is not a row to act on.

It renders every row, so it is for lists up to a few hundred items. A longer list is a `VirtualList`
(`new VirtualList({keyOf, render})`, then `setItems(items)`; its `selectedIndex` is the current row) beside a detail of
your own, in a `Splitter`.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (the height a view would give) against {oldLines} with
`ui.splitH` (third tab), which builds a row per compound by hand, has no keyboard and needs CSS of its own for the
selected row. u2 handles the current row and `aria-selected`, the arrow keys, the sash and its keyboard, and keeps
each detail it built.
