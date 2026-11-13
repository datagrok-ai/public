---
title: A panel toolbar and breadcrumbs
group: Navigation
order: 20
summary: Buttons, a toggle and a menu in one row that folds what does not fit into "…"; a path to climb back up.
live: yes
height: 200
---

`Toolbar` is one compact row for a panel's commands — buttons, toggles bound to a signal, menus — that folds
whole items into a "…" menu (More actions) when it runs out of width (narrow the window). It is one tab stop with
←/→ between the items, the toolbar pattern screen readers expect. `Breadcrumbs` shows a path; a click on a segment
hands you its index. A view's main controls go to its ribbon instead (see Inside Datagrok).

## What it takes

{lines} lines. u2 handles the overflow menu (re-measured on resize), roving focus, the toggle's `aria-pressed`,
the menus, and the breadcrumbs' current-page semantics.
