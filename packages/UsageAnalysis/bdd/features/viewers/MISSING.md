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

## Rules shared by every viewer (`viewer-rules-edge-cases.feature`)

The rows of the Charts viewers (sunburst, tree, word cloud) are in the Charts copy of the rules,
`packages/Charts/bdd/features/viewers/viewer-rules-edge-cases.feature`.

One scenario outline per rule, one Examples row per viewer; a row that breaks its rule is a known
failure in an Examples block of its own. Run on localhost 2026-10-29 (core master 575f6ae603): 57
tests, the known failures failing at their pinned steps.

### Steps requested

- **Min above Max is swapped** (bar chart and box plot value axis, line chart X and Y and the
  control limits, radar Min/Max map, word cloud text size and rotation; not the timelines, where an event
  whose End is before its Start is drawn as a point: `packages/Charts/bdd/features/viewers/timelines/`). No outline: no viewer reports the drawn range of its value axis. The
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
  logged; the context panel shows the missing name. A row of the removal outline, `@ticket-pending`.
- A layout naming a removed column puts that name back into the viewer's setting (scatter plot
  Color, PC plot columns, tree map split, sunburst and tree hierarchies); the viewers draw without
  the column and report no error. Not a defect (2026-10-30): the layout outline claims the name is
  kept. The forms viewer throws a TypeError on every render instead: public#4452, its own block.
- A viewer that picks another column after its column is removed (the box plot's category) is
  correct: the removal outline claims only that the removed column leaves the setting.
- **The tile viewer shows an empty value as a lane with no name** (lane names "x, y, "), not
  `(empty)`: core#1217.
- **A viewer's click filter stays after its column is removed** (bar chart, tree, sunburst): the
  removal does not refilter the table, and the rows come back on the next refilter (any filter
  change). Added to core#1174 (2026-10-30); pinned on it.
