# Tooltips — what is still missing

The features translated from TestTrack `Tooltips/` (seven md cases and their specs, six features:
uniform-default-tooltip is claimed by default-tooltip-visibility, whose first scenario hovers the same three
viewers with the same setup). What could not be said is listed here. Run on the local master stand,
2026-10-07.

## Resolved

- **The grid's own default tooltip shows nothing — by design** (Olesia, 2026-10-07). Since GROK-19901 (core
  7bf24ebcf5, the tooltip harmonization, changed it from "inherit from table") the grid's Show Tooltip defaults
  to "show custom tooltip" with an empty Row Tooltip, and the grid then shows no row tooltip, not even for the
  columns pushed out of sight. The TestTrack cases (`default-tooltip`,
  `default-tooltip-visibility`, `edit-tooltip`, `uniform-default-tooltip`, `tooltip-properties`) now switch the
  grid to "inherit from table" first, and grid-visible-columns-tooltip claims the default.
- **Tooltip > Hide switches off every tooltip of the table — by design** (Dmytro Kovalyov in the team thread,
  2026-10-07). The histogram's bin and the bar chart's bar tooltips go too, and so do the viewers' custom
  tooltips (`_isTooltipShownForViewer` in `tooltip.dart` reads the table's `.tooltip-visibility`). The cases
  `default-tooltip-visibility` (step 4) and `tooltip-properties` (step 11) were updated;
  default-tooltip-visibility and tooltip-properties claim it.
- **A custom tooltip with an empty Row Tooltip lists only the viewer's own data columns — by design**
  (Olesia, 2026-10-07; `tooltip.dart` has built it from `getExpandedColumns([])` since 2023): the scatter
  plot its X and Y (Data Values = Merge), the box plot nothing, not the table's tooltip columns.
  `tooltip-properties` step 7 was updated; tooltip-properties claims both.
- **A point of the line chart to hover** (2026-10-07 review). The line chart's status reports its first 20
  drawn points as `point of row <n>` (`line_chart_status.dart`; rows from 1, of the aggregated frame when
  aggregated); line-chart-aggregated-tooltip hovers the first one and claims the two aggregations in its
  tooltip (#2571).
- **A menu item that is not listed** (2026-10-07 review). `the open menu should not list` read the group's
  items the moment it hovered it, so a submenu not open yet read as an empty list and the negative passed;
  it now waits for the group to list something (mutation-tested: a listed item claimed absent fails).
- **The column selector of a form that rebuilds its rows** (2026-10-07 review). Edit Aggregated Tooltip redraws
  every row when one is added or changed, and a press on a selector at that moment opened nothing (1 run in
  10). The gesture presses again when no column grid opened.

## Missing steps and signals

- **The order of the tooltip's columns — a step, and a question.** The cases uniform-default-tooltip and
  edit-tooltip ask for "the same columns in the same order" on every viewer. `the tooltip should show columns {string}`
  uppercases, de-duplicates and sorts what it reads (`src/runtime/viewer-legend.ts`, `tooltipColumns`), so
  neither the order nor a column listed twice is claimed. Wanted: `Then the tooltip should show columns in
  the order {string}` — the first cells of `.d4-row-tooltip-table` rows as they come. Likely a product
  question too: with Data Values = Merge the scatter plot puts its axis columns first
  (`tooltip.dart`, `getExpandedColumns`), so its order differs from the box plot's and the grid's.
- **The trellis plot's own tooltip after Hide.** Its inner viewers' marks are not areas of the trellis, so
  default-tooltip-visibility leaves it out. The line chart is left out there too, for another reason: on
  demog-1000 it aggregates, and an aggregated line chart shows no tooltip until Edit Aggregated Tooltip
  configures one.
- **The choices of a property in the context panel.** `"Show Tooltip" property in context panel should offer
  the choice …` read an empty list once in three runs: the Dart property grid's choice editor puts its
  `<select>` in the cell only after the cell is clicked. tooltip-properties reads the default and picks the
  other two values instead. Wanted: `{element} should offer the choice` to open a lazy choice editor before reading it.
- **A menu item enabled or disabled.** No step claims it (`aria-disabled` on the item).
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
