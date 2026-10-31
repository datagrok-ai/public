@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot zoom and filter synchronization
  Zoom and Filter in all four modes: a zoom filters the table and Reset View gives it back, pack
  and zoom by filter does not leave the table stuck filtered, an external filter drives the
  viewport under zoom by filter, the filter panel's own reset undoes the plot's contribution (the
  status bar's "Filtered:" count appears with the zoom and goes with the reset), and
  Filter Out Invalid drops exactly the rows a logarithmic axis cannot draw; then the md's fifth and
  seventh scenarios — the axis type switched to logarithmic from the context menu (the property row
  is greyed out) on a datetime axis keeps the SEX card's 553 rows, and a large jitter on a
  logarithmic Y axis filters no row (the md's spgi-100 identifiers become WEIGHT by HEIGHT, whose
  values are all positive). Under filter by zoom
  the table's filter follows the viewport once the zoom animation has landed, so a claim on the
  filter is read after the zoom. Closed issues last: a window set with the range slider under
  Zoom and Filter = no action survives the filter panel opening, filtering and its reset
  (public#552); a plot shrunk to 40 by 40 adds no zoom filter of its own (public#1348); a
  logarithmic axis on an all-empty column under filter by zoom filters no row (public#1744). One journey on demog-1000, X = WEIGHT, Y = HEIGHT (872 of the
  1000 rows have a HEIGHT); every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows
    And all rows should pass the filter

  Scenario: A zoom filters the table and Reset View gives it back
    Then "Zoom and Filter" property of scatter plot viewer should be "filter by zoom"
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 1000 rows should pass the filter
    And scatter plot viewer should show fewer rows than before
    And the "rows selected" reading of scatter plot viewer should be the same as before
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: Pack and zoom by filter does not leave the table filtered
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 1000 rows should pass the filter
    When user sets "Zoom and Filter" property of scatter plot viewer to "pack and zoom by filter"
    And user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    When user sets "Zoom and Filter" property of scatter plot viewer to "filter by zoom"
    And user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: Under zoom by filter an external filter narrows the viewport
    When user sets "Zoom and Filter" property of scatter plot viewer to "zoom by filter"
    And user remembers the value range of scatter plot viewer
    And user remembers the "x axis span" reading of scatter plot viewer
    And user filters rows where "WEIGHT" is between 90.96 and 115.64
    Then 212 rows should pass the filter
    And the "x axis span" reading of scatter plot viewer should be lower than before
    And scatter plot viewer should show fewer rows than before
    And the "rows selected" reading of scatter plot viewer should be the same as before
    When user resets the filter
    Then all rows should pass the filter
    And the "x axis span" reading of scatter plot viewer should be as remembered
    And scatter plot viewer should show the remembered value range
    When user sets "Zoom and Filter" property of scatter plot viewer to "filter by zoom"
    Then no errors should have been logged

  Scenario: The filter panel's reset undoes what the zoom filtered
    When user opens the filter panel
    Then all rows should pass the filter
    And status bar should not show the text "Filtered:"
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 1000 rows should pass the filter
    And scatter plot viewer should show fewer rows than before
    And the "rows selected" reading of scatter plot viewer should be the same as before
    And status bar should show the text "Filtered:"
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And status bar should not show the text "Filtered:"
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A logarithmic axis alone filters nothing; Filter Out Invalid drops the non-positive rows
    When user adds a calculated column "AGE_SHIFT" with formula "${AGE} - 40"
    And user sets properties of scatter plot viewer:
      | X | WEIGHT    |
      | Y | AGE_SHIFT |
    Then scatter plot viewer should show 1000 rows
    And "Filter Out Invalid" property of scatter plot viewer should be "false"
    When user sets "Y Axis Type" property of scatter plot viewer to "logarithmic"
    Then all rows should pass the filter
    When user sets "Filter Out Invalid" property of scatter plot viewer to "true"
    Then 635 rows should pass the filter
    And the filter should pass exactly the rows where "AGE" is between 40.5 and 100
    When user sets "Filter Out Invalid" property of scatter plot viewer to "false"
    Then all rows should pass the filter
    When user sets properties of scatter plot viewer:
      | Y Axis Type        | linear |
      | Filter Out Invalid | true   |
    Then all rows should pass the filter
    And scatter plot viewer should show 1000 rows
    When user sets properties of scatter plot viewer:
      | Filter Out Invalid | false  |
      | X                  | WEIGHT |
      | Y                  | HEIGHT |
    And user removes "AGE_SHIFT" column
    And user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: The axis type switched from the menu on a datetime axis leaves the filter as it was
    When user picks "STARTED" in the "x" column selector of scatter plot viewer
    And user adds a categorical filter on "SEX" keeping "F"
    Then 553 rows should pass the filter
    And "X Axis Type" property of scatter plot viewer should be "linear"
    When user picks "Properties... > X Axis > X Axis Type > Logarithmic" from the context menu of scatter plot viewer
    Then scatter plot viewer should have repainted
    And "X Axis Type" property of scatter plot viewer should be "logarithmic"
    And 553 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    When user picks "Properties... > X Axis > X Axis Type > Linear" from the context menu of scatter plot viewer
    Then "X Axis Type" property of scatter plot viewer should be "linear"
    When user adds a categorical filter on "SEX" keeping "F, M"
    And user sets "X" property of scatter plot viewer to "WEIGHT"
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: A large jitter on a logarithmic axis filters no row
    When user sets properties of scatter plot viewer:
      | jitterSize  | 16 |
      | jitterSizeY | 17 |
    Then scatter plot viewer should have repainted
    When user sets "Y Axis Type" property of scatter plot viewer to "logarithmic"
    Then scatter plot viewer should have repainted
    And the "y axis min" reading of scatter plot viewer should differ from before
    And all rows should pass the filter
    And scatter plot viewer should show 872 rows
    When user sets properties of scatter plot viewer:
      | Y Axis Type | linear |
      | jitterSize  | 0      |
      | jitterSizeY | 0      |
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

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
    Then 359 rows should pass the filter
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
