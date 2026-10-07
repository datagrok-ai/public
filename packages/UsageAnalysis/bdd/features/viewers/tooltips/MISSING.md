# Tooltips — what is still missing

The features translated from TestTrack `Tooltips/` (seven md cases and their specs). Written with today's
steps only (no library, binding or core change); what could not be said is listed here. Run on the local
master stand, 2026-10-07.

## Resolved

- **The grid's own default tooltip shows nothing — by design** (Olesia, 2026-10-07). Since GROK-19901 the
  grid's Show Tooltip defaults to "show custom tooltip" with an empty Row Tooltip, and the grid then shows no
  row tooltip, not even for the columns pushed out of sight. The TestTrack cases (`default-tooltip`,
  `default-tooltip-visibility`, `edit-tooltip`, `uniform-default-tooltip`, `tooltip-properties`) now switch the
  grid to "inherit from table" first, and grid-visible-columns-tooltip claims the default.
- **Tooltip > Hide switches off every tooltip of the table — by design** (Dmytro Kovalyov in the team thread,
  2026-10-07). The histogram's bin and the bar chart's bar tooltips go too, and so do the viewers' custom
  tooltips (`_isTooltipShownForViewer` in `tooltip.dart` reads the table's `.tooltip-visibility`). The cases
  `default-tooltip-visibility` (step 4) and `tooltip-properties` (step 11) were updated;
  default-tooltip-visibility and tooltip-properties claim it.
- **A custom tooltip with an empty Row Tooltip lists only the viewer's own data columns — by design**
  (Olesia, 2026-10-07): the scatter plot its X and Y (Data Values = Merge), the box plot nothing, not the
  table's tooltip columns. `tooltip-properties` step 7 was updated; tooltip-properties claims both.

## Missing steps and signals

- **The order of the tooltip's columns — a step, and a question.** The cases uniform-default-tooltip and
  edit-tooltip ask for "the same columns in the same order" on every viewer. `the tooltip should show columns {string}`
  uppercases, de-duplicates and sorts what it reads (`src/runtime/viewer-legend.ts`, `tooltipColumns`), so
  neither the order nor a column listed twice is claimed. Wanted: `Then the tooltip should show columns in
  the order {string}` — the first cells of `.d4-row-tooltip-table` rows as they come. Likely a product
  question too: with Data Values = Merge the scatter plot puts its axis columns first
  (`tooltip.dart`, `getExpandedColumns`), so its order differs from the box plot's and the grid's.
- **A point of the line chart to hover.** line-chart-aggregated-tooltip's own subject, the aggregated tooltip
  over a dot of a split chart (#2571), needs a place a pointer can rest on: the line chart's status reports
  `view`, `plot`, `chart N` and the axes, no point, and the chart's centre raises nothing. Wanted in
  `line_chart_status.dart`: `point <n> of <series>` (or `marker of row <r>` as the scatter plot has), the
  marker's box from the last frame. Then:

  ```gherkin
    When user hovers over the "point 1 of chart 1" area of line chart viewer
    Then tooltip should contain text "concat unique(Stereo Category)"
    And tooltip should contain text "min(Average Mass)"
  ```
- **The trellis plot's and the line chart's own tooltips after Hide.** default-tooltip-visibility has them in
  the view; neither reports a place that raises a tooltip (the trellis inner viewers' points are not areas of
  the trellis). Same need as above.
- **The choices of a property in the context panel.** `"Show Tooltip" property in context panel should offer
  the choice …` read an empty list once in three runs: the Dart property grid's choice editor puts its
  `<select>` in the cell only after the cell is clicked. tooltip-properties reads the default and picks the
  other two values instead. Wanted: `{element} should offer the choice` to open a lazy choice editor before reading it.
- **A menu item that is not listed.** `the open menu should not list "Tooltip > Hide"` re-hovers the group and
  reads its items at once (`viewer-menus.ts`, `menuLabels` with no wanted item), so a submenu not open yet
  reads as an empty list and the negative passes. The features pair it with a positive on the same group
  (`… should list "Tooltip > Show Custom"`) just before. Wanted: the negative waits for the group to show at
  least one item. Likewise no step claims a menu item enabled (`aria-disabled`).
- **The Edit Tooltip dialog's list has three columns** (type, name, box). The list steps read the names and the
  boxes only; a reading of the list's own grid columns (`column order` of the grid inside a dialog) would
  claim it. Low value.
- **The Add aggregation button of Edit Aggregated Tooltip has no label.** Its Dart name is cut
  (`button-Add-aggregation-to-calculate-and-show-in-to`), so the feature reaches it as `first button`. A
  caption or an aria-label on the button in core would let it be named.

## Not translated, by choice

- **default-tooltip step 4, widening the context panel** to push the last column out of sight — it moves
  nothing the widened column does not; grid-visible-columns-tooltip pushes the columns off by the column
  resizer.
- **line-chart-aggregated-tooltip's "in-viewer filter"** — named in the case's title only; none of its steps
  sets one.

## Not translated, by rule

- **What the tooltip's renderer draws** (molecules in SPGI's tooltip) — a renderer, not a feature
  (`libraries/bdd/CLAUDE.md`). edit-tooltip and the others use demog-1000 for this reason.
- **The `.tooltip` and `.tooltip-visibility` tags the old specs read** — echoes of what the gesture wrote;
  the features claim the tooltip a hover shows instead.
