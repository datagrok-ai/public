# Layout edge cases — what the features cannot say yet

What `layout-edge-cases.feature` leaves out, written for the edge-case features with no library, binding or code
change.

## Still open

**Questions for owners**

1. **A layout of a grid alone: suggested for every table, or only for tables with the grid's columns?**
   `help/visualize/view-layout.md` ("Layout suggestions") says a grid-only layout has no columns of interest and is
   suggested for any table. `TableView.getViewerColumns` (`table_view.dart`) says that with no other viewer the
   grid's columns are used. On localhost (2026-10-26) a grid-only layout saved on a demog table is not suggested
   for cars. One of the two is to change before a scenario can claim it.
2. **A layout renamed to an empty name in the context panel**: refused, or put back?
3. **A layout applied (by a click on its card) to a table that lacks one of its columns**: the viewer shows the
   missing-column placeholder as a dashboard does, or the card is not offered at all? The suggestion rule keeps such
   a card out of the pane, but the gallery (View > Layout > Open Gallery) lists every layout.

**API tests to write**

1. **`ApiTests/src/dapi/layouts.ts`, "get applicable" (`:225`) and "filter" (`:238`)** assert
   `layouts.length >= 0`, which always holds. "get applicable" should find the layout the test saved for its
   table, and not find a layout whose viewer uses a column the table lacks; "filter" should find the saved layout
   by its name (its `[0]` at `:236` throws on an empty list, so it fails only by accident).

## Won't do

- Rename followed by Edit tags (core#689) and the panel title after a rename (core#688): written as known failures
  in `layout-names-tags.feature` on branch `bdd/layout-names-tags`; not repeated here.
