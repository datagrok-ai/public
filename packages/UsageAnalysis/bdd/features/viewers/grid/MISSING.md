# Grid — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026, against the grid's
"Automation surface" (`core/client/d4/lib/src/viewers/grid/CLAUDE.md`), the TestTrack cases
`files/TestTrack/Viewers/Grid/grid-ui.md`, `grid.md`, `grid-rows-select-filter-navigate.md` and the
library's rules (`public/libraries/bdd/CLAUDE.md`). Added in that round: `grid-pin-if-editable.feature`,
`grid-extracted-rows-project.feature` (GROK-19717), and the exact row count in `grid-columns.feature`
(326 after the popup's min handle reaches bin 10). Control+F and the block drag across the cells stay
out of scope, as `grid-rows.feature` records.

## Still open

- **The header histogram and the stats rows** (`grid-ui.md` "Column Stats — Visual Verification").
  The grid reports `stats rows` (which rows are drawn), not what they show: needs a
  `text of <stat> of <col>` reading (the string the stats row paints, `grid_special_rows.dart`) to
  claim "the Min row shows the minimum AGE". Same for the header histogram (`showColumnHistograms`):
  nothing says it was drawn for a column; a `header histogram <col>` area would.
- **Form designer field positions** (`grid-ui.md` "Summary Columns — Form Designer"): a field box
  dragged in the designer "settles where it was dropped" and comes back on Edit. The designer's field
  boxes carry no names and report no positions (`tile-viewer.ts` reads only which fields exist).
  Needs a name per field box (`form-field-<column>`) and its position in the designer's status.
- **The hamburger popup's categorical filter** (`grid-ui.md` "Context Panel — Column Hamburger Menu
  Inline Filter", steps 4-5): `grid-columns.feature` drives the histogram filter of the popup; for a
  string column the popup holds a categorical filter widget whose rows have no hit areas of their own
  (the filter panel's `category <label> of <caption>` areas are the panel's, not the popup's). Needs
  the popup's filter to publish `category <label>` areas like the panel's card does.

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed grid issues into `grid-regressions.feature` (localhost, master).

### Regressions found 26 Oct 2026

- **A cell whose text is cut off shows its whole value in the tooltip (public#989, public#996).**
  Measured twice on localhost: on demog-1000 the USUBJID column dragged 50 px narrower wraps
  `X0273T21000300003` and cuts its second line off; hovering `cell 1 of USUBJID` shows a tooltip box
  with no text in it (`tooltip should be visible` passes, `tooltip should contain text
  "X0273T21000300003"` fails: no visible tooltip holds the text). An uncut cell (`cell 1 of RACE`)
  shows no tooltip at all. A numeric column does not cut off: HEIGHT narrowed by 25 px drops a digit
  (`174.7`) and shows the same empty box. `grid_core.showCellTooltip` passes the cut-off text as
  `preText`, and the "show custom tooltip" branch of `tooltip.dart` (the grid's default since
  GROK-19901) does not append it. Expected (the issue): the whole value in the tooltip.

  ```gherkin
  @realizes:public#989 @realizes:public#996
  Scenario: A cell whose text is cut off shows its whole value in the tooltip (public#989)
    Given user is logged in
    And user opens demog-1000 dataset
    Then the "text of cell 1 of USUBJID" reading of grid should be "X0273T21000300003"
    When user drags the "column resizer USUBJID" area of grid by 50 pixels to the left
    Then the "column width of USUBJID" reading of grid should be lower than before
    When user hovers over the "cell 1 of USUBJID" area of grid
    Then tooltip should be visible
    And tooltip should contain text "X0273T21000300003"
    When user moves the pointer away from grid
    Then no errors should have been logged
  ```

### Steps requested

- **Toolbox search, then Extract Rows on the rows it found (public#3963).** The table view's search box
  is not reachable by a phrase: with the toolbox shown and its Search pane expanded, `"…" input in
  toolbox` finds no input ("inputs in toolbox: none on the page", probed 26 Oct 2026). Needs a step
  `user searches {string} in the toolbox` (types into the table view's search box and presses Enter),
  then:

  ```gherkin
  @realizes:public#3963
  Scenario: The rows the toolbox search finds are the rows Extract Rows extracts (public#3963)
    Given user is logged in
    And user opens demog-1000 dataset
    And the toolbox pane is shown
    When user searches "RACE = Asian" in the toolbox
    Then the filter should pass exactly the rows where "RACE" is "Asian"
    And 15 rows should pass the filter
    When user picks "Select > Extract Rows" from the top menu
    Then the table should have 15 rows
    And the table should have no rows where "RACE" is "Caucasian"
    And no errors should have been logged
  ```
  (15 Asian rows on demog-1000; the menu path of Extract Rows to be checked when the step exists.)
