---
title: Settings per item, beside the list
group: Panels & layout
order: 32
summary: Three calculations with a check each, the current one's settings beside them, out of reach until it is checked.
live: yes
height: 260
---

`MasterDetail` is a short list beside the detail of its current item, split by a sash the user drags. A row has
a label, an optional `check` (two-way: a calculation runs or not) and a `hint` that follows what the detail
holds. The detail is built the first time its item is current and kept, so what was typed in it survives a move
to another item; `disabled` keeps it in view, dimmed and out of reach, with the reason over it. ↑/↓, Home and
End move the current item, Space checks it. It fills its host: the sample's CSS gives the host a height. A long
list of records to sort and group is a `CatalogTable` instead; thousands of rows beside one record's details are
a `VirtualList` in a `Splitter` (A list with details).

## What it takes

{lines} lines and {cssLines} of CSS (the host's height). u2 handles the rows, the current one and its keyboard,
`aria-selected` and `aria-checked`, the sash, the kept details and the inert, dimmed detail with its reason.
