# Color coding editor — what the features cannot say yet

Written during the color scheme editor round of 14 Oct 2026 (TestTrack `Viewers/color-coding-ui.md`
steps 9-12, the help's "Color code columns", and the closed Jira and GitHub tickets on color coding).
Translation changed no code: each item is a missing step, name or signal with the scenario it would
let us write.

## 1. A category's colour in the categorical editor (public#3413, GROK-10437)

Each category of "Color-coding: <col>" is a `.d4-color-code-item` (reachable as `"M" item`), but the
picker opens on its colour bar, `.d4-legend-item-bar`, which has no name; a click on the item itself
opens nothing. Needed: a name on the bar (e.g. `name="color-bar-<category>"`) or a `bar` part on the
`item` kind for these lists. With it, `user picks the color "<hex>" in the color picker` (which
clicks the picker's `color-<hex>` swatch) finishes the scenario:

```gherkin
  Scenario: A colour picked for a category in the editor is the colour the grid paints it
    When user picks "Color Coding > Categorical" from the context menu of the "header SEX" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header SEX" area of grid
    And user clicks on the colour bar of "M" item in "Color-coding: SEX" dialog
    And user picks the color "<a swatch hex>" in the color picker
    Then the categorical color of "M" in "SEX" column should be "<hex>"
    And the "color of cell 4 of SEX" reading of grid should be "#<hex>"
    When user selects "Off" in Type input in "Color-coding: SEX" dialog
    Then "SEX" column should have no color coding
```

The same bar is what GROK-10437 (the editor showed black for categories set by hex) would read: the
bar's colour against the column's categorical colour.

## 2. The conditional rows and the stops of "Edit color scheme"

The conditional editor's rows (a `MapInput`: a pattern box with the placeholder "(no value)" and a
colour cell) and the stops of "Edit color scheme" (a value and a colour per row) have no names; the
features reach only the rows' Add row / Remove row buttons by position. Not reachable today:
- a pattern typed into a row, an invalid pattern's `d4-invalid` mark and its "Invalid pattern"
  tooltip;
- a bound value that colours the cells equal to it (public#1967);
- Below min / Above max: they are `div-Below-min` / `div-Above-max`, a name no kind matches, so the
  "only one picker is open" case (GROK-17738) cannot be driven;
- a stop's colour changed in the scheme.

Needed: names on the rows' inputs (`input-host-pattern-<n>` / `input-host-color-<n>`), and a kind for
the two `div-…` colour inputs.

## 3. The text colour of a cell

`Apply to = text` is claimed by the column's text flag (`the text of X column should be color-coded`
reads `isTextColorCoded`), and `color of cell <r> of <c>` is the renderer's back colour, which is the
coded colour in text mode too. So nothing checks that the renderer paints the letters and leaves the
cell's background alone (public#2784). Needed: a `text color of cell <r> of <c>` reading from the
`TextCellRenderer`'s text/back colour split.

Nor does a step claim the flag unset (`the text of X column should not be color-coded`). Apply to =
background does not remove the `.color-coding-text` tag: `GridColumn.isTextColorCoded = false` sets it
to null (`grid_column.dart`, `setTag` without `removeIfNull`), which the column facts read as the
string `"null"`, so `"X" column should have tag ".color-coding-text" equal to ""` fails while the flag
is correctly off. Claiming `"null"` would pin that detail, so the features switch back to background
without a claim, and the layout scenario cannot claim the flag gone before the layout brings it back.

## 4. Answered

- **No Cancel in the editors.** "Edit color scheme" has only OK (`cancelText: ''`,
  `color_coding.dart`) and the conditional editor only CLOSE; every change applies at once. Olesia
  Pavlenko, 14 Oct 2026: this is the intended behaviour, and `color-coding-ui.md` step 10 (Cancel
  reverts the scheme) and GROK-8444 (Cancel in the conditional editor reverts) describe a dialog that
  no longer exists, so the TestTrack case is stale. The features claim the current behaviour instead:
  a removed stop repaints the column before OK (`color-coding-editor-linear.feature`).
- **Off keeps the text flag.** `Color Coding > Off` leaves the column's text flag set, so a later
  Linear comes back on the text: intended (Olesia Pavlenko, 14 Oct 2026). The features switch Apply
  to back to background before Off.

## 5. Dropped by rule

- What a layout or a project stores for a coloring, the absolute-scale tags before and after a
  layout: ApiTests `Grid: Color Coding` (`colorCoding.api`).
- What a cell renderer draws, a correlation plot's mirrored cells (GROK-20919): pixel checks of a
  renderer.
- The legend's own picker, its Cancel and a second pick (github-3132): `legend/legend-across-viewers.feature`;
  a conditional coloring's legend: `guides/scatter-plot-conditional-colors.feature`.
