# Layout edge cases — what the features cannot say yet

What `layout-edge-cases.feature` leaves out, written for the edge-case features with no library, binding or code
change.

## Still open

**Decided and written**

The answers of 2026-10-27 are scenarios of `layout-edge-cases`: a layout of a grid alone is suggested only for a
table that has the grid's columns (`TableView.getViewerColumns`; fixed in core#1179 on 2026-10-28, a known failure until the
stand has it: the lookup read `dataFrame[name]` instead of `dataFrame[colName]`); a layout renamed to an empty name keeps its name. View > Layout > Open Gallery opens the same Layout
suggestions panel, under the same rule, so a layout missing a column of the table is offered nowhere.

**API tests to write**

1. **`ApiTests/src/dapi/layouts.ts`, "get applicable" (`:225`) and "filter" (`:238`)** assert
   `layouts.length >= 0`, which always holds. "get applicable" should find the layout the test saved for its
   table, and not find a layout whose viewer uses a column the table lacks; "filter" should find the saved layout
   by its name (its `[0]` at `:236` throws on an empty list, so it fails only by accident).

## Won't do

- Rename followed by Edit tags (core#689) and the panel title after a rename (core#688): written as known failures
  in `layout-names-tags.feature` on branch `bdd/layout-names-tags`; not repeated here.
