@journey @viewers @realizes:viewers.heat-map
Feature: Heat map scroll sliders, rectangle selection and the Table property
  The gestures the old checklist kept as manual. A drag of a scroll slider's handle narrows that
  axis only: `y scroll span` (or `x scroll span`) drops below 1 while the other stays at 1, and
  the map shows it — taller rows for the vertical one, fewer columns and the headers of the
  columns left of the window gone for the horizontal one; a double-click on the slider puts its
  axis back. Navigation's Alt-drag zoom covers both axes at once (heat-map-navigation.feature).
  A Shift-drag over the AGE band selects the rows the rectangle spans in the table — the inner
  80% of the band's height, rows 100 to 892 of the thousand — and Escape clears the selection.
  Bound to another table through the Table property, the map draws that table's rows and
  columns (cars, 30 rows, 17 columns, each row tall enough for per-cell areas), and bound back it
  draws demog-1000 again.

  Background:
    Given user is logged in
    And user opens cars dataset
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "x scroll span" reading of heat map viewer should be 1
    And the "y scroll span" reading of heat map viewer should be 1
    And heat map viewer should be painted

  Scenario: Dragging the vertical slider's handle narrows the rows, and a double-click resets them
    When user drags the "y scroll max handle" area of heat map viewer by 300 pixels to the up
    Then the "y scroll span" reading of heat map viewer should be between 0.5 and 0.9
    And the "x scroll span" reading of heat map viewer should be 1
    And the "row height" reading of heat map viewer should be higher than before
    And the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should have repainted
    When user double-clicks on the "y scroll slider" area of heat map viewer
    Then the "y scroll span" reading of heat map viewer should be 1
    And the "row height" reading of heat map viewer should be lower than before
    And no errors should have been logged

  Scenario: Dragging the horizontal slider's handle narrows the columns, and a double-click resets them
    When user drags the "x scroll min handle" area of heat map viewer by 200 pixels to the right
    Then the "x scroll span" reading of heat map viewer should be between 0.4 and 0.9
    And the "y scroll span" reading of heat map viewer should be 1
    And the "columns shown" reading of heat map viewer should be lower than before
    And heat map viewer should not have a "header USUBJID" area
    And heat map viewer should have a "header SEVERITY" area
    And heat map viewer should have repainted
    When user double-clicks on the "x scroll slider" area of heat map viewer
    Then the "x scroll span" reading of heat map viewer should be 1
    And the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should have a "header USUBJID" area
    And no errors should have been logged

  Scenario: A Shift-drag selects the rows it spans and Escape clears them
    Given user clears the row selection
    When user drags a selection box over the "column AGE" area of heat map viewer
    Then rows 100 to 892 should be selected
    And heat map viewer should show a selection highlight
    When user presses Escape
    Then no rows should be selected
    And no errors should have been logged

  Scenario: The Table property binds the map to another table and back
    When user sets "table" property of heat map viewer to "cars"
    Then heat map viewer should be bound to table "cars"
    And the "rows" reading of heat map viewer should be 30
    And the "columns shown" reading of heat map viewer should be 17
    And heat map viewer should have a "cell 30 of model" area
    And the "text of cell 1 of model" reading of heat map viewer should be "alfaromeo"
    And heat map viewer should not have a "column AGE" area
    And heat map viewer should have repainted
    When user sets "table" property of heat map viewer to "demog-1000"
    Then heat map viewer should be bound to table "demog-1000"
    And the "rows" reading of heat map viewer should be 1000
    And the "columns shown" reading of heat map viewer should be 11
    And the "column order" reading of heat map viewer should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And the "row height" reading of heat map viewer should be between 0 and 8
    And heat map viewer should have a "column AGE" area
    And no errors should have been logged
