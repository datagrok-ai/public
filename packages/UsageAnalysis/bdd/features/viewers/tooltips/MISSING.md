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
- **The order of the tooltip's columns** (2026-10-12): `the tooltip should show columns in the order {string}`
  reads the rows as they come. The scatter plot's axis-first order under Data Values = Merge is by design (the
  lead, 2026-10-11): `edit-tooltip` claims each viewer's own order — the scatter plot HEIGHT, WEIGHT, AGE, SEX,
  the box plot and the grid AGE, SEX, WEIGHT (mutation-checked: a swapped order fails).
- **The choices of a property in the context panel** (2026-10-12): `{element} should offer …` opens a property
  grid row's lazy choice editor before reading it, so the one-in-three empty read is gone.
- **A menu item enabled or disabled** (2026-10-12): `the open menu should list {string} as enabled` / `as disabled`
  read the item's `aria-disabled`.
- **The Add aggregation button of Edit Aggregated Tooltip** (2026-10-12): the button of an empty list is named
  and labelled "Add aggregation" in core (`aggregated_columns_input.dart`); `line-chart-aggregated-tooltip` names
  it instead of `first button`.

## Still open

Nothing planned.

## Won't do

- **The trellis plot's own tooltip after Hide.** Its inner viewers' marks are not areas of the trellis, so
  default-tooltip-visibility leaves it out; the line chart too, which aggregates on demog-1000 and shows no
  tooltip until Edit Aggregated Tooltip configures one.
- **The Edit Tooltip dialog's list has three columns** (type, name, box): a reading of a dialog's inner grid
  columns for one low-value line.

### Not translated, by choice

- **default-tooltip step 4, widening the context panel** to push the last column out of sight — it moves
  nothing the widened column does not; grid-visible-columns-tooltip pushes the columns off by the column
  resizer.
- **line-chart-aggregated-tooltip's "in-viewer filter"** — named in the case's title only; none of its steps
  sets one.

### Not translated, by rule

- **What the tooltip's renderer draws** (molecules in SPGI's tooltip) — a renderer, not a feature
  (`libraries/bdd/CLAUDE.md`). edit-tooltip and the others use demog-1000 for this reason.
- **The `.tooltip` and `.tooltip-visibility` tags the old specs read** — echoes of what the gesture wrote;
  the features claim the tooltip a hover shows instead.

## Edge cases (`tooltips-edge-cases.feature`)

Written for `tooltips-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data shapes,
columns removed or renamed under the open viewer, layouts naming a missing column, property bounds,
closing the viewer. Nothing here asks for a library, binding or code change made by that feature; it
lists what the feature could not claim and why.

### Edge cases: Still open

#### 1. Steps requested

None.

#### 2. Decided

1. **A tooltip column renamed** is shown under its new name (the rename rule of
   `viewer-rules-edge-cases.feature`). The table tag keeps column names (`tooltip.dart:574-577`) and
   nothing renames them, so the column drops out: core#1215, a known failure in
   `tooltips-edge-cases.feature`.
2. **Ten columns ticked, seven shown**: the cut stays, and the tooltip says that part of the
   columns was dropped (the silent-cut rule). Today an explicit selection is cut to `colAmount` = 7
   (`tooltip.dart:590`) with no sign of it. Steps requested: the note has no wording yet, and no step
   reads the lines of a tooltip; then: tick ten columns, hover a row, claim seven value lines and the
   note.

3. **Edit Aggregated Tooltip with an added row** (Add aggregation, OK): the row takes a column by
   default and is saved as `avg(AGE)`: right as it is (2026-10-30), claimed in
   `tooltips-edge-cases.feature` on a demog-1000 line chart (X = AGE).

#### 3. Questions for owners

None.

#### 3. API or Dart tests to write

None.

### Edge cases: Won't do

Nothing.
