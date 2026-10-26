@viewers @realizes:viewers.grid
Feature: Grid edge cases
  Guards for the grid at the edges of its settings and its data, on demog-1000 (AGE in row 4 is 45,
  SEVERITY in row 3 is None and in row 4 Medium). An HTML summary column keeps rendering after a
  column it shows is removed, and the grid scrolls to its last row with no error (TestTrack
  `Viewers/Grid/grid.md`, "Summary Column Surviving a Removed Source Column", the GROK-19942 contract).
  Two scenarios describe defects: Unpin Column with Frozen Columns above the column count
  (core#1137), and text typed into a numeric cell reported as a successful edit, which has no ticket
  yet (`@ticket-pending`).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows

  @edge-cases
  Scenario: An HTML summary column keeps rendering after a column it shows is removed
    When user picks "Add > Summary Columns > Default HTML Form" from the context menu of the "cell 2 of USUBJID" area of grid
    And user clicks on All link in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "cell type of html" reading of grid should be "html"
    When user picks "Remove" from the context menu of the "header CONTROL" area of grid
    Then the table should not have a column "CONTROL"
    When user clicks on the "cell 2 of AGE" area of grid
    And user presses Control+End
    Then grid should have a "cell 1000 of html" area
    And the "cell type of html" reading of grid should be "html"
    And no errors should have been logged

  @edge-cases @known-failure:37 @realizes:core#1137
  Scenario: Unpin Column with Frozen Columns above the column count logs no error
    core#1137: with Frozen Columns 20 on an 11-column table every column counts as frozen, and Unpin
    Column inserts the column at position 19 of a shorter list, which throws "RangeError: Value not in
    range". Remove the tag when core#1137 is fixed and this scenario passes.
    When user sets "Frozen Columns" property of grid to "20"
    And user picks "Pin > Unpin Column" from the context menu of the "header AGE" area of grid
    Then no errors should have been logged

  @edge-cases @ticket-pending @known-failure:51
  Scenario: Text typed into a numeric cell is not reported as an edit
    The defect: the grid opens a text editor on the integer AGE column, the typed "abc" does not parse
    and the value stays, yet the grid fires its cell-value-edited event as for an accepted edit (and
    shows the user nothing). Remove the tag when the rejected text no longer fires the event.
    Given user listens for "d4-grid-cell-value-edited" event on grid
    When user double-clicks on the "cell 4 of AGE" area of grid
    And user presses Control+A in cell editor
    And user types "abc" into cell editor
    And user presses Enter
    Then cell editor should be hidden
    And the value of "AGE" column in row 4 should be "45"
    And "d4-grid-cell-value-edited" event should not have fired on grid
