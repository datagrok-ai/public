# Grid — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

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
