@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot regressions — filtering, zoom, selection and labels
  Closed scatter plot issues turned into claims on what the plot draws: labels for selected rows
  the filter took out are not drawn (public#3189), a viewer Filter survives a new calculated column
  (public#1156), a zoom set with the range slider under Zoom and Filter = no action survives the filter panel opening, filtering
  and its reset (public#552; a wheel zoom is animated past the plot's settle, so the slider sets the window), a
  plot shrunk to a corner does not filter rows by its tiny viewport (public#1348), Mouse Drag =
  Select turns a plain drag into a selection and Pan back into a pan (public#2089), Show Selected
  Rows = false draws no selection highlight (public#3613), the labels of the selected rows all come
  back after a zoom in and out (public#3176), a click on a marker after a zoom makes its row current
  (public#801), a selection box takes the markers whose Color value is empty (public#1174), the range
  slider moves a logarithmic axis (public#376), and a logarithmic axis on an all-empty column under
  filter by zoom filters no row (public#1744). One journey on demog-1000, X = WEIGHT, Y = HEIGHT
  (872 of the 1000 rows have a HEIGHT); every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows
    And no rows should be selected

  # public#3189
  @realizes:public#3189
  Scenario: Labels for selected rows are not drawn for the rows the filter took out (public#3189)
    When user sets properties of scatter plot viewer:
      | Label Columns   | AGE      |
      | Show Labels For | Selected |
    And user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    When user selects rows where "SEX" is "M"
    Then the "labels shown" reading of scatter plot viewer should be at least 1
    When user selects rows where "SEX" is "F"
    Then 553 rows should be selected
    And 0 selected rows should pass the filter
    And the "labels shown" reading of scatter plot viewer should be 0
    When user clears the row selection
    And user adds a categorical filter on "SEX" keeping "F, M"
    And user sets properties of scatter plot viewer:
      | Label Columns   |     |
      | Show Labels For | All |
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#1156
  @realizes:public#1156
  Scenario: A viewer Filter keeps its rows when a calculated column is added (public#1156)
    When user sets "Filter" property of scatter plot viewer to '${RACE} in ["Asian", "Black"]'
    Then scatter plot viewer should show 38 rows
    When user adds a calculated column "AGE2" with formula "${AGE} * 2"
    Then the table should have a column "AGE2"
    And scatter plot viewer should show 38 rows
    And scatter plot viewer should report no error
    When user removes "AGE2" column
    And user sets "Filter" property of scatter plot viewer to ""
    Then scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#552
  @realizes:public#552
  Scenario: A zoom under Zoom and Filter = no action survives filtering in the filter panel (public#552)
    When user sets properties of scatter plot viewer:
      | Zoom and Filter    | no action |
      | axesFollowFilter   | false     |
    And user drags the min handle of the "x" range slider of scatter plot viewer by 30 pixels
    Then scatter plot viewer should show a narrower value range than before
    When user remembers the value range of scatter plot viewer
    And user opens the filter panel
    Then scatter plot viewer should show the remembered value range
    When user adds a range filter on "AGE" from 20 to 40
    Then fewer than 1000 rows should pass the filter
    And scatter plot viewer should show the remembered value range
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And scatter plot viewer should show the remembered value range
    When user sets properties of scatter plot viewer:
      | Zoom and Filter    | filter by zoom |
      | axesFollowFilter   | true           |
    And user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#1348
  @realizes:public#1348
  Scenario: A plot shrunk almost to nothing filters no row, and the panel's reset gives every row back (public#1348)
    When user opens the filter panel
    And user resizes scatter plot viewer to 40 by 40
    And user adds a range filter on "WEIGHT" from 100 to 170
    Then 143 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    When user restores the size of scatter plot viewer
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#2089
  @realizes:public#2089
  Scenario: Mouse Drag = Select makes a plain drag select, and Pan makes it pan again (public#2089)
    When user sets "Zoom and Filter" property of scatter plot viewer to "no action"
    And user remembers the "x axis min" reading of scatter plot viewer
    And user sets "mouseDrag" property of scatter plot viewer to "Select"
    And user drags across the "view" area of scatter plot viewer
    Then some rows should be selected
    And the "x axis min" reading of scatter plot viewer should be as remembered
    When user selects no rows
    And user sets "mouseDrag" property of scatter plot viewer to "Pan"
    And user drags across the "view" area of scatter plot viewer
    Then the "x axis min" reading of scatter plot viewer should not be as remembered
    And no rows should be selected
    When user picks "Reset View" from the context menu of scatter plot viewer
    And user sets "Zoom and Filter" property of scatter plot viewer to "filter by zoom"
    Then all rows should pass the filter
    And the "x axis min" reading of scatter plot viewer should be as remembered
    And no errors should have been logged

  # public#3613
  @realizes:public#3613
  Scenario: Show Selected Rows = false draws the selected markers without the highlight (public#3613)
    When user selects rows where "SEX" is "F"
    Then scatter plot viewer should show a selection highlight
    When user sets "showSelectedRows" property of scatter plot viewer to "false"
    Then 553 rows should be selected
    And scatter plot viewer should show no selection highlight
    When user sets "showSelectedRows" property of scatter plot viewer to "true"
    Then scatter plot viewer should show a selection highlight
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  # public#3176
  @realizes:public#3176
  Scenario: Every label of the selected rows comes back after a zoom in and out (public#3176)
    When user sets properties of scatter plot viewer:
      | Zoom and Filter | no action |
      | Label Columns   | AGE       |
      | Show Labels For | Selected  |
    And user selects the first 5 rows
    Then the "labels shown" reading of scatter plot viewer should be 5
    When user remembers the value range of scatter plot viewer
    And user scrolls the mouse wheel up 4 times over the "view" area of scatter plot viewer
    Then scatter plot viewer should show a narrower value range than before
    When user scrolls the mouse wheel down 4 times over the "view" area of scatter plot viewer
    Then the "labels shown" reading of scatter plot viewer should be 5
    When user clears the row selection
    And user picks "Reset View" from the context menu of scatter plot viewer
    And user sets properties of scatter plot viewer:
      | Label Columns   |                |
      | Show Labels For | All            |
      | Zoom and Filter | filter by zoom |
    Then the "labels shown" reading of scatter plot viewer should be 0
    And all rows should pass the filter
    And no errors should have been logged

  # public#801
  @realizes:public#801
  Scenario: After a zoom a click on a marker makes its row current (public#801)
    When user scrolls the mouse wheel up 2 times over the "view" area of scatter plot viewer
    Then scatter plot viewer should show fewer rows than before
    When user clicks on the first "marker of row" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be the current row
    And the "hovered row" reading of scatter plot viewer should not be "1"
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#1174
  @realizes:public#1174
  Scenario: A selection box takes the markers whose Color value is empty (public#1174)
    When user adds a calculated column "RACE_F" with formula "if(${SEX} == \"F\", ${RACE}, null)"
    And user sets "Color" property of scatter plot viewer to "RACE_F"
    And user drags a selection box over the "view" area of scatter plot viewer
    Then some rows where "SEX" is "M" should be selected
    And some rows where "SEX" is "F" should be selected
    When user clears the row selection
    And user sets "Color" property of scatter plot viewer to ""
    And user removes "RACE_F" column
    Then no rows should be selected
    And no errors should have been logged

  # public#376
  @realizes:public#376
  Scenario: The range slider moves a logarithmic axis (public#376)
    When user sets properties of scatter plot viewer:
      | Zoom and Filter | no action   |
      | X Axis Type     | logarithmic |
    And user remembers the "x axis min" reading of scatter plot viewer
    And user drags the min handle of the "x" range slider of scatter plot viewer by 30 pixels
    Then the "x axis min" reading of scatter plot viewer should be higher than remembered
    And the "x axis min" reading of scatter plot viewer should be a finite number
    And the "x axis max" reading of scatter plot viewer should be a finite number
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then the "x axis min" reading of scatter plot viewer should be as remembered
    When user sets properties of scatter plot viewer:
      | X Axis Type     | linear         |
      | Zoom and Filter | filter by zoom |
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#1744
  @realizes:public#1744
  Scenario: A logarithmic axis on an all-empty column under filter by zoom filters no row (public#1744)
    When user adds a calculated column "NO_VALUES" with formula "if(${AGE} > 1000, ${WEIGHT}, null)"
    Then "NO_VALUES" column should have missing values
    And "Zoom and Filter" property of scatter plot viewer should be "filter by zoom"
    When user sets properties of scatter plot viewer:
      | Y Axis Type | logarithmic |
      | Y           | NO_VALUES   |
    Then scatter plot viewer should show 0 rows
    And all rows should pass the filter
    When user sets properties of scatter plot viewer:
      | Y           | HEIGHT |
      | Y Axis Type | linear |
    And user removes "NO_VALUES" column
    And user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged
