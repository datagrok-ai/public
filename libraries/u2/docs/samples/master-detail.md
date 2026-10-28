---
title: A list with details
group: Panels & layout
order: 30
summary: Two thousand compounds in a virtual list beside the selected one's details, split by a sash.
live: yes
old: ui.splitH today
---

The master-detail shape: a `VirtualList` (only the visible rows are in the DOM, so the count does not matter)
beside a details pane, in a `Splitter` the user drags or moves with the arrow keys. The selection is a signal
(`selectedIndex`), the details are built once and follow it, and `show()` swaps in the empty state while nothing
is selected. Click a row, or focus the list and use ↑/↓. A detail pane on its own is a `tableFromMap`: muted keys,
values that wrap, a value that is an element (here a live `span`) shown as it is, and no hover tint — it is not a
row to act on.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (the heights a view would give) against {oldLines} with
`ui.splitH` (third tab) — about the same length, but it builds a DOM row per compound, has no keyboard and needs
CSS of its own for the selected row. u2 handles the virtualization, the list's keyboard and `aria-selected`, the
sash and its keyboard, and the empty state.
