# Tree map — what the features cannot say yet

Written during the heat map / 3D scatter plot / tree map gap round of 19 Oct 2026, against
the viewer's "Automation surface" (`core/client/d4/lib/src/viewers/tree_map/CLAUDE.md`) and the library's
rules (`public/libraries/bdd/CLAUDE.md`). `tree-map-properties.feature` holds what could be said:
Show Selected Rows, Row Source Selected, a bar chart filtering the map (public#935), Auto Layout and the
selectors (public#2131), and the categorical gate.

## Still open

- Picking a column in an on-viewer split selector: the selectors are plain `<select>`s with no name
  (`TreeMapCore._addSplit`); needs `split selector <n>` names, or the area-based pick learning a native `<select>`.
- `{widget} should report the error {string}` on the tree map: its status reports the gate as the `error` reading only;
  the status should fill `error:` too.
- No selection highlight on the RACE colours: the category colours fall in the selection hue; needs a `bands drawn`
  reading.

## Won't do

- What a person sees on a table with no categorical column: Add Viewer does not offer the tree map (the lead,
  2026-10-19); the feature claims the gate an API-added map shows.
- The selection band's height: the features claim the band present or absent, paired with a repaint.
- What a layout stores for Show Selected Rows / Auto Layout: an API test's; `tree-map-color-size-selection.feature`
  claims the one the md asks for.
- Color or Size removed (core#1135), NaN in a Size or Color column, a renamed split column, no rows, 1000 categories:
  xamgle `d4 | Viewers | Column removal | Tree map…` and `d4 | Viewers | Non-finite values | Tree map`, the rules every
  viewer shares (`../MISSING.md`).
