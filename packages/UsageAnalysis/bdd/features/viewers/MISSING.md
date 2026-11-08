# Cross-viewer features — what cannot be said yet

The features in this folder rather than in a viewer's own: `viewers-docking.feature`,
`viewers-layout-download.feature` (View > Layout > Download gives a `.layout` file; what it holds is the layout API's
test) and `molecule-categories.feature`; and where the rules every viewer shares are checked.

## Rules shared by every viewer

Each rule is tested once where it lives; a viewer's feature claims it only where its own gesture is the subject.

- A value past a property's declared min or max is brought to the bound (core#1200): prop_gen `property_test.dart`
  "Property.clamp keeps a number inside min and max…", xamgle `d4 | Viewers | Fixed issues | A property keeps a value
  inside its declared bounds` (the panel's input, `setOptions`, a layout), ApiTests `AI: Viewers: JsViewer` "a property
  keeps a value inside its declared min and max".
- No rows, every row deleted, every column removed, empty, constant and non-finite values: xamgle `d4 | Viewers | Data
  Manipulation | <viewer>` and `d4 | Viewers | Non-finite values | <viewer>` (every d4 viewer); Charts `Viewers: no rows`.
- A removed column leaves the settings that named it: `DataFrameViewer` (`look.auto`), xamgle `d4 | Viewers | Column
  removal | …`; Charts `subscribeToColumnRemoval` and `Viewers`.
- A renamed column is followed by every setting that named it, the Filter formula included (core#1213): ApiTests
  `AI: Viewers: Shared rules`.
- A layout naming a missing column is applied and the viewer drawn: xamgle `d4 | Viewers | Invalid settings | Trellis
  plot reports a split column the table lacks…` and `d4 | Viewers | Fixed issues | Scatter plot reports its status after
  a layout with missing columns`, ApiTests `Shell: ViewLayout` (public#2867).
- A click filter goes when its viewer is closed: a gesture, claimed in each filtering viewer's own feature
  (`bar-chart-interaction`, `pie-chart-interaction`, `pc-plot-filter`, `calendar-filter`,
  `filter-panel/compose-with-viewers` for the scatter plot's zoom, `sunburst`, `tree-on-click-filter`,
  `trellis-plot-click-actions`).

## Resolved

- core#1174: a click filter goes when the clicked column is removed — the sunburst and the tree (Charts; claimed in
  `sunburst.feature` and `tree-on-click-filter.feature`), the bar chart and the filter-panel cards (d4; their Dart
  tests, see `bar-chart/MISSING.md` and `filter-panel/MISSING.md`).

## Still open

- public#4452 — the forms viewer throws a TypeError on every render after a layout naming a removed column.
- core#1239 — the markup viewer renders HTML from a layout or a project as trusted (Mode Html, Markup off); what may stay
  is the owner's call.
- Min above Max is swapped (bar chart and box plot value axis, line chart, radar, word cloud): no viewer reports the
  window it drew; needs `value axis min` / `value axis max` readings.
- The splitter between two docked viewers: dock_spawn's splitters carry no name and report no box; needs a name per
  splitter, then `… should be wider than before` after a drag.
- Dropping the downloaded `.layout` back onto the view: needs a step that drops the file just downloaded.
- The range sliders over molecule categories: the box plot reports no slider area, the trellis plot no `first x
  category` reading.

## Won't do

- Every viewer property set through the Context Panel instead of the API step: the property grid's editors are generic,
  claimed once per kind.
- `(empty)` names of empty categories per viewer (Sankey nodes, pivot keys, tree map leaves): one category formatter;
  the sunburst and the tree claim it.
- A viewer announcing the part it cuts (pivot groups, Sankey categories, sunburst levels, tooltip columns): the caps are
  by design.
- core#714's dock drop preview box: cosmetic, for a closed bug; it needs a dock_spawn name for one claim.
- The per-viewer rule scenarios of the 2026-11-08 round (out of bounds, no rows, deleted rows, removed or renamed
  columns, layouts naming a missing column): the tests above.
