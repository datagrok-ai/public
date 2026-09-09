@journey @viewers @realizes:viewers.histogram
Feature: Histogram bin selection, row markers and mouse-over
  A click on a bin selects exactly the rows whose value falls in it — both edges included — Control
  adds a bin and a Shift drag spans the bins it crosses; the selected share of each bin is painted
  over the bar and disappears under a selected row source. The current row and the mouse-over row
  are dots on the baseline, the second one driven from another viewer, and hovering a bin shows the
  bin's tooltip and dims every other bin behind the rows of the group. One journey on demog-1000
  with a histogram of AGE; bin 8 is [42.85, 46.4] — AGE 43 to 46, 99 rows.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value | AGE |
    And user resizes histogram viewer to 500 by 400
    Then histogram viewer should show 1000 rows
    And histogram viewer should have a "bin 8" area
    And histogram viewer should show no selection highlight
    And histogram viewer should not have a "selected bin 8" area

  Scenario: A click on a bin selects the rows of that bin
    Given user listens for "d4-histogram-select-bins" event on histogram viewer
    When user clicks on the "bin 8" area of histogram viewer
    Then "d4-histogram-select-bins" event should have fired on histogram viewer
    And only rows where "AGE" is one of "43, 44, 45, 46" should be selected
    And 99 rows should be selected
    And histogram viewer should show a selection highlight
    And histogram viewer should have a "selected bin 8" area
    And the "selected bin 8" area of histogram viewer should contain the color "#FF8C00"
    And histogram viewer should not have a "selected bin 1" area
    And no errors should have been logged

  Scenario: Control adds a bin to the selection
    When user clicks on the "bin 9" area of histogram viewer holding Control
    Then only rows where "AGE" is one of "43, 44, 45, 46, 47, 48, 49" should be selected
    And 180 rows should be selected
    And histogram viewer should show more selection highlight than before
    And histogram viewer should have a "selected bin 9" area
    And no errors should have been logged

  Scenario: Clearing the selection takes the overlay with it
    When user clears the row selection
    Then no rows should be selected
    And histogram viewer should show no selection highlight
    And histogram viewer should not have a "selected bin 8" area
    And histogram viewer should not have a "selected bin 9" area
    And no errors should have been logged

  Scenario: A Shift drag selects every bin it crosses
    When user drags a selection box from the "bin 6" area to the "bin 7" area of histogram viewer
    Then only rows where "AGE" is one of "36, 37, 38, 39, 40, 41, 42" should be selected
    And histogram viewer should show a selection highlight
    And histogram viewer should have a "selected bin 6" area
    And histogram viewer should have a "selected bin 7" area
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: A selected row source shows the selection instead of overlaying it
    When user selects rows where "RACE" is "Asian"
    And user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    And histogram viewer should show 553 rows
    And histogram viewer should show a selection highlight
    When user sets "Row Source" property of histogram viewer to "Selected"
    Then histogram viewer should show 15 rows
    And histogram viewer should show no selection highlight
    And only rows where "RACE" is "Asian" should be selected
    When user sets "Row Source" property of histogram viewer to "All"
    Then histogram viewer should show 1000 rows
    And histogram viewer should show a selection highlight
    When user sets "Row Source" property of histogram viewer to "Filtered"
    And user resets the filter
    And user clears the row selection
    Then all rows should pass the filter
    And no rows should be selected
    And no errors should have been logged

  Scenario: The current row is a dot on the baseline
    When user makes row 8 current
    Then row 8 should be current
    And histogram viewer should have a "current row marker" area
    And the "current row marker" area of histogram viewer should contain the color "#38B738"
    When user sets "Show Current Row" property of histogram viewer to "false"
    Then histogram viewer should not have a "current row marker" area
    And histogram viewer should have repainted
    When user sets "Show Current Row" property of histogram viewer to "true"
    Then histogram viewer should have a "current row marker" area
    And no errors should have been logged

  Scenario: The mouse-over row marker follows the grid
    Given histogram viewer should not have a "mouse over row marker" area
    When user hovers over the "cell 8 of AGE" area of grid
    Then histogram viewer should have a "mouse over row marker" area
    And the "mouse over row marker" area of histogram viewer should contain the color "#AAAAAA"
    When user moves the pointer away from grid
    Then histogram viewer should not have a "mouse over row marker" area
    And no errors should have been logged

  Scenario: Hovering a bin shows its tooltip and dims the other bins
    Given user listens for "d4-histogram-mouse-over-bins" event on histogram viewer
    When user hovers over the "bin 8" area of histogram viewer
    Then "d4-histogram-mouse-over-bins" event should have fired on histogram viewer
    And the "bin 1" area of histogram viewer should have repainted
    And histogram viewer should have repainted by at least 500 pixels
    And tooltip should contain text "total: 99"
    And tooltip should contain text "filtered: 99"
    When user moves the pointer away from histogram viewer
    Then no errors should have been logged

  Scenario: A mouse-over row group in another viewer repaints the bins
    When user adds a bar chart viewer with:
      | Split | RACE |
    And user takes a snapshot of histogram viewer
    And user hovers over the "bar Asian" area of bar chart viewer
    Then histogram viewer should have repainted by at least 500 pixels
    When user moves the pointer away from bar chart viewer
    And user moves the pointer away from histogram viewer
    And user sets "Show Mouse Over Row Group" property of histogram viewer to "false"
    And user hovers over the "bar Black" area of bar chart viewer
    Then histogram viewer should not have repainted
    When user moves the pointer away from bar chart viewer
    And user sets "Show Mouse Over Row Group" property of histogram viewer to "true"
    And user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent
    And no errors should have been logged
