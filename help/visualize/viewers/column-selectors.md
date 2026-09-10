---
title: "Column Selector"
description: Use the combo box control for choosing and previewing a column to plot on a viewer axis.
keywords:
  - column picker
  - combo box
  - axis selector
  - column preview
  - formula
---

A combo box for selecting a column

**Mouse-over** to see column summary, as well as a viewer preview.
**Drag a column** into a selector to change it.
**Start typing** to filter columns.
**Esc** to reset column filter, or hide popup.
**Enter** to accept current choice.
**Space** to pop up column selector.

## Using a formula instead of a column

Any selector that drives a viewer — an axis, color, size, markers, whiskers, split — can be given a
formula rather than an existing column. Hover the viewer and click the **plus** icon next to the
selector to open [Add New Column](../../transform/add-new-column.md), then write an expression such as
`${AGE} * 2`.

The result is stored as a hidden column, so it does not clutter the grid, and it is shared: once
defined, the same formula can be picked from any other selector on that table. To change it, select it
and click the **pencil** icon — the formula is edited in place, so every viewer using it updates. To
keep the original and add a variant, select a regular column first and use the plus icon again.

Formula columns are listed at the end of the column picker.

See also:

* [Viewers](../viewers/viewers.md)
* [Table View](../table-view-1.md)
