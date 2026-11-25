@platform
Feature: Selecting rows by their row headers
  A run of rows selected the way a person selects it in the grid: a click on the first row's header makes that
  row the selection and the current row, a Shift+click on the last one's adds every row between them. Filmed,
  the two presses are two stops.

  Scenario: A click and a Shift+click on the row headers select the rows between them
    Given user is logged in
    And user opens cars dataset
    When user selects rows 3 to 7 by their row headers
    Then rows 3 to 7 should be selected
