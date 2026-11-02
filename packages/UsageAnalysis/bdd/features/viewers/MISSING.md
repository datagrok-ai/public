# Cross-viewer features — what cannot be said yet

Written during the viewer gaps round of 23 Oct 2026 for the features that sit in
this folder rather than in a viewer's own: `viewers-docking.feature`,
`viewers-layout-download.feature` (added in that round: View > Layout > Download gives a `.layout`
file; what it holds is the layout API's test) and `molecule-categories.feature`.

## Still open

- **Properties set through the API step, not the Context Panel** (every viewer's TestTrack case: "In Context
  Panel > X set Y"). Most viewer features change a property with `user sets {string} property of {widget}`,
  which writes the look directly; the property-grid editor the user actually clicks is exercised only in a
  few places (`user selects ... in ... property in context panel`, `user checks ... property in context
  panel`). Not rewritten in this round: switching every feature to the property-grid gestures is a sweep over
  others' features. The gestures exist; what is missing is the decision to use them by default.
- **The splitter between two docked viewers** (`viewers-docking-ui.md`: "drag the splitter between
  them"). dock-spawn's splitter bars carry no name and report no box, so neither a drag nor "the
  left viewer got wider" can be aimed. Needs a name per splitter (`dock-splitter-<n>`, or named by
  the two panels it divides) in `core/client/libs/dock_spawn`; the claim is then
  `histogram viewer should be wider than before` from a `user drags {element} by N and 0 pixels`.
- **Dropping a downloaded `.layout` back onto the view** (`viewers-docking-ui.md`, last step): a
  file dragged in from the desktop. `user drops the {string} file of the project onto {element}`
  takes a fixture file, not the file just downloaded; a step that drops the downloaded file would
  close the round trip.
- **The range sliders over molecule categories on the box plot and the trellis plot**
  (`rendering-structures-on-the-axes-ui.md`, #2330). The bar chart and the scatter plot are walked
  in `molecule-categories.feature`. The box plot reports no slider area (`box_plot/CLAUDE.md`), so
  its category slider cannot be aimed at. The trellis plot reports its `x scroll handle`, and a drag
  of it moves the categories on screen, but every reading stays the same (`x categories shown` 5,
  `x scroll handle share` 0.05 on spgi-100 by Structure) and the label areas are named by the
  molblock text — no phrase can spell them. Needs a `first x category` reading (its row index or
  position) on the trellis, and `x range slider` / handle areas on the box plot:

  ```gherkin
  Scenario: Dragging the trellis plot's category scroll handle over the molecules moves the window
    When user adds a trellis plot viewer with:
      | xColumnNames | Structure |
    Then the "x scroll handle share" reading of trellis plot viewer should be 0.05
    When user drags the "x scroll handle" area of trellis plot viewer by 100 pixels to the right
    Then the "first x category" reading of trellis plot viewer should be higher than before
    And the "blank cells" reading of trellis plot viewer should be 0
    And trellis plot viewer should report no error
  ```

## Rules shared by every viewer

A filter left by a viewer: the rows a viewer's click filtered out come back when the viewer is
closed, or when the column of the clicked category is removed. Out of bounds: a value typed in the
context panel outside the property's declared bounds is not taken; the property keeps its value.
The panel's number input marks such a value invalid (`number_input.dart`, `addMinMaxValidators`),
but the input writes every change into the property (`InputBase.bind`).

The checks of these rules for the sunburst, the tree and the word cloud are steps of their
journeys: `sunburst/sunburst.feature`, `tree/tree-on-click-filter.feature`, `tree/tree.feature` and
`word-cloud/word-cloud-appearance.feature`. For the bar chart, calendar, histogram, pie chart and
trellis plot: `bar-chart/bar-chart-tooltip.feature`, `bar-chart/bar-chart-selection-overlays.feature`,
`bar-chart/bar-chart.feature`, `calendar/calendar.feature`, `histogram/histogram-range-filter.feature`,
`pie-chart/pie-chart.feature` and `trellis-plot/trellis-plot-pick-up-color-and-scroll.feature`.

A viewer that breaks a rule carries the check as the known failure of its ticket. Run on localhost
2026-10-29 (core master 575f6ae603): the known failures fail at their pinned steps.

### Steps requested

- **Min above Max is swapped** (bar chart and box plot value axis, line chart X and Y and the
  control limits, radar Min/Max map, word cloud text size and rotation; not the timelines, where an event
  whose End is before its Start is drawn as a point, requested from ApiTests in `timelines/MISSING.md`). No outline: no viewer reports the drawn range of its value axis. The
  line chart's `x axis min` / `x axis max` are its range slider's bounds, and the box plot, the bar
  chart, the radar and the word cloud report no axis range at all. Wanted: `value axis min` /
  `value axis max` readings of the drawn window, then rows that set Min 60 and Max 20 on AGE and
  claim the window between 20 and 60.
- **Empty categories as `(empty)`** on the Sankey, the pivot table's key cells and the tree map: no
  reading lists their node, key or leaf names.
- **No current row highlights nothing** on the radar: no reading of the highlighted row.

### Security, for the owner

- **The markup viewer renders HTML from a layout or a project as trusted.** With Mode Html (or Auto
  over content starting with `<`) and Markup off, `markup_viewer_core.dart` puts the content into the
  DOM with `NodeTreeSanitizer.trusted`, so `<script>` and `<img onerror=...>` from someone else's
  layout or project run in the viewer's session. The d4 rule is `NodeTreeSanitizer.trusted` only for
  constants. No scenario: whether the content should be sanitized, and what may stay (styles,
  links), is the owner's call: core#1239.

### Found on the run, no ticket yet

- **The bar chart's Stack keeps the name of a removed column.** Stacking is dropped and nothing is
  logged; the context panel shows the missing name. No ticket yet.
- A layout naming a removed column puts that name back into the viewer's setting (scatter plot
  Color, PC plot columns, tree map split, sunburst and tree hierarchies); the viewers draw without
  the column and report no error. Not a defect (2026-10-30): the name is kept. The forms viewer
  throws a TypeError on every render instead: public#4452, still open.
- A viewer that picks another column after its column is removed (the box plot's category) is
  correct: the check is only that the removed column leaves the setting.
- **The tile viewer shows an empty value as a lane with no name** (lane names "x, y, "), not
  `(empty)`: core#1217.
- **A viewer's click filter stays after its column is removed** (bar chart, tree, sunburst): the
  removal does not refilter the table, and the rows come back on the next refilter (any filter
  change). Added to core#1174 (2026-10-30); pinned on it.

### API or Dart tests to write

demog-1000, readings and settings of the viewer as named; every check also asserts that nothing is
logged and the viewer reports no error.

- **No rows**: with every row filtered out (AGE between 500 and 600) the viewer shows 0 rows and
  draws nothing; after the reset it shows 1000 rows and draws again. Bar chart (split RACE), box
  plot (value AGE), correlation plot, heat map, histogram (value AGE), matrix plot, network diagram
  (node 1 RACE), PC plot (AGE, HEIGHT, WEIGHT; lines drawn), pie chart (RACE; slices), pivot table
  (aggregated rows), scatter plot (X AGE), tile viewer (lanes RACE; tiles), tree map (RACE, SEX;
  leaves), trellis plot (Y RACE).
- **Every row deleted**: what the viewer draws goes from at least 1 to 0. Bar chart (bars), box
  plot, calendar, heat map, line chart (X AGE; markers drawn), matrix plot, network diagram, PC
  plot, pie chart, pivot table.
- **A column removed** leaves the setting that named it: scatter plot Color RACE and Size HEIGHT,
  tile viewer lanes RACE, tree map Color AGE and Size WEIGHT, box plot Category 1 SEX.
- **A column renamed** is followed by the setting: scatter plot X AGE → AGE_YEARS and Color RACE →
  ORIGIN, pie chart category RACE → ORIGIN, PC plot columns HEIGHT → TALL (AGE, TALL, WEIGHT), tile
  viewer lanes RACE → ORIGIN, tree map split SEX → GENDER (RACE, GENDER), trellis plot Y RACE → ORIGIN.
- **A Filter formula renamed**: a scatter plot Filter `${AGE} > 44` reads `${YEARS} > 44` after AGE
  is renamed to YEARS. core#1213, still open: the rename handler renames column settings only.
- **A layout naming a removed column** is applied, the viewer is shown and the setting keeps the
  name: scatter plot Color RACE, PC plot columns AGE, HEIGHT, WEIGHT without HEIGHT, tree map split
  RACE, SEX without SEX.
- **Out of bounds through the API**: histogram Bins 201 and 0 keep 20, pie chart Max Radius 600
  keeps 150, trellis plot Tiles Per Row 0 keeps 4. core#1200, still open: the value is written.
- **An empty value** is the `(empty)` category: a tile viewer with lanes on a column holding an
  empty value shows 3 tiles and a lane named `(empty)`. core#1217, still open: the lane has no name.

Moved out of BDD (no user gesture); requested from ApiTests.
