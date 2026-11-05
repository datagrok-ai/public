---
title: A catalog of records
group: Lists & tables
order: 20
summary: Campaigns to sort, group, choose columns of and act on — every choice remembered.
live: yes
height: 420
old: ui.table today
---

`CatalogTable` is the list of records a person comes back to — campaigns, profiles, templates. Describe the
columns (`value` sorts and groups, `render` draws), the groupings, the row actions and where to keep the user's
choices; it owns the order, the groups and the shown columns. Double-click a header to sort (a date column sorts
as dates, codes in natural order: KHL-2 before KHL-12); the icons beside the heading choose the columns and the
grouping; a group folds on a click. A hovered row takes the platform grid's hover tint and shows its icon actions;
a right click opens all of them. An action's `run` returns a promise or nothing.

## What it takes

{lines} lines against {oldLines} with `ui.table` (second tab) — which sorts on a header click (codes in text
order: KHL-12 before KHL-2) by rebuilding the table, and has no groups, no column or group menus, no context
menu and no memory. u2 handles the sorting (`aria-sort`, Enter on a focused header), the groups and their
counts, both menus, the row actions and the context menu, the pinned name column in a narrow pane, the footer,
and saving it all under one key.
