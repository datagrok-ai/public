@journey @viewers @realizes:viewers.grid
Feature: Grid Tags column
  What a Tags column draws: a tag in the rows whose source column carries the flag, and nothing
  over the rest of the grid - the cells of the other columns keep their own background and the
  other summary columns keep their own drawing. The claims are written as the product should
  behave and the scenario carries `@known-failure` for GROK-20888: on dev a Tags column paints the
  grid over in one colour, so the scenario is expected to fail and its passing will say the ticket
  is closed and the tag has to go. It lives apart from `summary-columns.feature` because the
  library inverts `@known-failure` only inside a `@journey`. A flood passes any "is painted"
  claim, so the claims are the white background of an unmarked Tags cell, two colours in a marked
  one, and the white background of a neighbouring column and of the Sparklines column; the values
  of the neighbouring columns are read as text beside them. One journey on demog-1000 (CONTROL is
  false in rows 1 and 2, true in row 3).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    Then grid should show 1000 rows

  @known-failure
  Scenario: A Tags column marks the rows that carry the flag and leaves the rest of the grid alone (GROK-20888)
    When user picks "Add > Summary Columns > Sparklines" from the context menu of the "cell 2 of USUBJID" area of grid
    And user picks "Add > Summary Columns > Tags" from the context menu of the "cell 2 of Sparklines" area of grid
    Then the "cell type of Tags" reading of grid should be "tags"
    And the "cell 1 of Tags" area of grid should contain the color "#FFFFFF"
    And the "cell 3 of Tags" area of grid should be painted in at least 2 colors
    And the "cell 1 of Sparklines" area of grid should contain the color "#FFFFFF"
    And the "cell 1 of SEVERITY" area of grid should contain the color "#FFFFFF"
    And the "text of cell 1 of SEVERITY" reading of grid should be "High"
    And the "text of cell 1 of CONTROL" reading of grid should be "false"
    And the "text of cell 3 of CONTROL" reading of grid should be "true"
    And no errors should have been logged
