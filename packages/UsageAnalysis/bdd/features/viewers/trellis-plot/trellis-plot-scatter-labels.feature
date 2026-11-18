@journey @viewers @realizes:viewers.trellis-plot @realizes:viewers.scatter-plot
Feature: Scatter plot labels in trellis cells
  A trellis starts its scatter plots on Show Labels For = FilteredSelected, and a scatter plot inside a cell
  labels only the selected rows it draws whatever Show Labels For says, at most 20 of them, with Display Labels Auto or Always; Never shows none. The cell opened full
  screen is a standalone scatter plot again: it labels by its Show Labels For, with no cap. A click on empty
  space in a cell, or in the cell opened full screen ("full screen view"), keeps the selection, which the same
  click on a standalone scatter plot clears. A cell's count is
  the trellis's "cell labels shown <cell>" reading, passed on from the scatter plot inside, and the full-screen
  one is "full screen labels shown". One journey on demog-1000 split by SEX and RACE (the scatter plot inside
  on HEIGHT and WEIGHT), labelled by USUBJID. Nothing is labelled until rows are selected. AGE 20-25 selects
  1 drawn row of F | Black, 3 of F | Other, none of M | Black and 28 of F | Caucasian, which shows 20; SEX F
  selects the 17 drawn rows of F | Black and the 429 of F | Caucasian (51 of its 480 rows have no HEIGHT or
  WEIGHT), which the cell caps at 20 and its full-screen view labels all. Each scenario clears the selection
  and puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    And user sets "labelColumnNames" inner property of trellis plot viewer to "USUBJID"
    Then the "cells" reading of trellis plot viewer should be 8
    And "showLabelsFor" inner property of trellis plot viewer should be "FilteredSelected"
    And the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 0
    And the "cell labels shown F | Other" reading of trellis plot viewer should be 0

  Scenario: A cell labels the selected rows it draws, whatever Show Labels For says
    When user selects rows where "AGE" is between 20 and 25
    Then the "cell labels shown F | Black" reading of trellis plot viewer should be 1
    And the "cell labels shown F | Other" reading of trellis plot viewer should be 3
    And the "cell labels shown M | Black" reading of trellis plot viewer should be 0
    And the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 20
    When user sets "showLabelsFor" inner property of trellis plot viewer to "Filtered"
    Then the "cell labels shown F | Other" reading of trellis plot viewer should be 3
    And the "cell labels shown M | Black" reading of trellis plot viewer should be 0
    When user sets "showLabelsFor" inner property of trellis plot viewer to "All"
    Then the "cell labels shown F | Other" reading of trellis plot viewer should be 3
    When user selects no rows
    Then the "cell labels shown F | Other" reading of trellis plot viewer should be 0
    When user sets "showLabelsFor" inner property of trellis plot viewer to "FilteredSelected"
    Then the "cell labels shown F | Black" reading of trellis plot viewer should be 0
    And no errors should have been logged

  Scenario: A cell shows at most 20 labels, with Display Labels Auto or Always, and none with Never
    When user selects rows where "SEX" is "F"
    Then the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 20
    And the "cell labels shown F | Black" reading of trellis plot viewer should be 17
    And the "cell labels shown M | Caucasian" reading of trellis plot viewer should be 0
    When user sets "displayLabels" inner property of trellis plot viewer to "Always"
    Then the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 20
    And the "cell labels shown F | Black" reading of trellis plot viewer should be 17
    When user sets "displayLabels" inner property of trellis plot viewer to "Never"
    Then the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 0
    And the "cell labels shown F | Black" reading of trellis plot viewer should be 0
    When user sets "displayLabels" inner property of trellis plot viewer to "Auto"
    Then the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 20
    When user selects no rows
    Then the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 0
    And no errors should have been logged

  Scenario: The cell opened full screen labels by its own Show Labels For, with no cap, and keeps the selection
    When user selects rows where "SEX" is "F"
    And user sets "showLabelsFor" inner property of trellis plot viewer to "Selected"
    Then the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 20
    When user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    And user clicks on the "full screen icon" area of trellis plot viewer
    Then "SEX: F, RACE: Caucasian" dialog should be visible
    And the "full screen labels shown" reading of trellis plot viewer should be 429
    When user listens for "d4-trellis-plot-inner-viewer-clicked" event on trellis plot viewer
    And user clicks on the "top left corner of full screen view" area of trellis plot viewer
    Then "d4-trellis-plot-inner-viewer-clicked" event should have fired on trellis plot viewer
    And only rows where "SEX" is "F" should be selected
    And the "full screen labels shown" reading of trellis plot viewer should be 429
    When user presses Escape
    Then "SEX: F, RACE: Caucasian" dialog should be absent
    And trellis plot viewer should not report a "full screen labels shown" reading
    When user selects no rows
    And user sets "showLabelsFor" inner property of trellis plot viewer to "FilteredSelected"
    Then the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 0
    And no errors should have been logged

  Scenario: A click on empty space in a cell keeps the selection a standalone scatter plot clears
    Given user listens for "d4-trellis-plot-inner-viewer-clicked" event on trellis plot viewer
    When user selects rows where "SEX" is "F"
    And user clicks on the "top left corner of cell body F | Caucasian" area of trellis plot viewer
    Then "d4-trellis-plot-inner-viewer-clicked" event should have fired on trellis plot viewer
    And only rows where "SEX" is "F" should be selected
    And the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 20
    When user adds a scatter plot viewer
    And user clicks on the "empty space" area of scatter plot viewer
    Then no rows should be selected
    And the "cell labels shown F | Caucasian" reading of trellis plot viewer should be 0
    When user clicks on close icon of scatter plot viewer
    Then no errors should have been logged
