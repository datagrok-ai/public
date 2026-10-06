@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot keyboard navigation
  The arrow keys walk the current cell through the grid, in the order the categories are drawn in —
  x from left to right, y from top to bottom — and a round trip of four arrows comes back where it
  started. Under On Click = Filter the arrow carries the filter to the cell it lands on, and Escape
  gives every row back. The keys go to the viewer's own root, which is focusable, so no click on
  the charts grid is needed to make them work. One journey on demog-1000 with SEX by RACE.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the cells of trellis plot viewer should be 2 wide and 4 tall

  Scenario: Four arrows walk the current cell around a square and back
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "F | Caucasian"
    When user presses ArrowRight in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "M | Caucasian"
    When user presses ArrowDown in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "M | Other"
    When user presses ArrowLeft in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "F | Other"
    When user presses ArrowUp in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "F | Caucasian"
    And no errors should have been logged

  Scenario: Every arrow announces the cell it moved to
    Given user listens for "d4-trellis-plot-current-cell-changed" event on trellis plot viewer
    When user clicks on the "cell F | Asian" area of trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "F | Asian"
    When user presses ArrowRight in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "M | Asian"
    When user presses ArrowDown in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "M | Black"
    When user presses ArrowLeft in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "F | Black"
    And "d4-trellis-plot-current-cell-changed" event should have fired on trellis plot viewer
    And no errors should have been logged

  Scenario: An arrow carries the filter to the cell it lands on
    When user sets "On Click" property of trellis plot viewer to "Filter"
    And user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 480 rows should pass the filter
    When user presses ArrowRight in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "M | Caucasian"
    And 416 rows should pass the filter
    And no rows where "SEX" is "F" should pass the filter
    And no errors should have been logged

  Scenario: Escape gives every row back
    When user presses Escape in trellis plot viewer
    Then all rows should pass the filter
    And the "current cell" reading of trellis plot viewer should be ""
    When user sets "On Click" property of trellis plot viewer to "None"
    Then "On Click" property of trellis plot viewer should be "None"
    And no errors should have been logged
