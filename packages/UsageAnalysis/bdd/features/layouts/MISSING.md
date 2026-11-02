# Layout edge cases — what the features cannot say yet

What the layout edge cases, now @edge-cases scenarios of `layout-names-tags.feature`, leave out, written for them with no library, binding or code
change.

## Still open

**Decided and written**

The answers of 2026-10-27 are @edge-cases scenarios of `layout-names-tags`: a layout of a grid alone is suggested only for a
table that has the grid's columns (`TableView.getViewerColumns`, fixed in core#1179: the lookup read
`dataFrame[name]` instead of `dataFrame[colName]`); a layout renamed to an empty name keeps its name. View > Layout > Open Gallery opens the same Layout
suggestions panel, under the same rule, so a layout missing a column of the table is offered nowhere.

**API tests to write**

1. **`ApiTests/src/dapi/layouts.ts`, "get applicable" (`:225`) and "filter" (`:238`)** assert
   `layouts.length >= 0`, which always holds. "get applicable" should find the layout the test saved for its
   table, and not find a layout whose viewer uses a column the table lacks; "filter" should find the saved layout
   by its name (its `[0]` at `:236` throws on an empty list, so it fails only by accident).

## Won't do

- Rename followed by Edit tags (core#689) and the panel title after a rename (core#688): claimed by the scenarios
  of `layout-names-tags.feature` ("A second layout saved and renamed takes a tag through Edit tags... on its card", and
  the title of the context panel after the rename); not repeated as edge cases.
