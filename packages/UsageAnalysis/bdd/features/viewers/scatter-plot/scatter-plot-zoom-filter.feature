@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot zoom and filter synchronization
  Zoom and Filter in all four modes: a zoom filters the table and Reset View gives it back, pack
  and zoom by filter does not leave the table stuck filtered, an external filter drives the
  viewport under zoom by filter, the filter panel's own reset undoes the plot's contribution, and
  Filter Out Invalid drops exactly the rows a logarithmic axis cannot draw. Under filter by zoom
  the table's filter follows the viewport once the zoom animation has landed, so a claim on the
  filter is read after the zoom. One journey on demog-1000, X = WEIGHT, Y = HEIGHT (872 of the
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
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 1000 rows should pass the filter
    And scatter plot viewer should show fewer rows than before
    And the "rows selected" reading of scatter plot viewer should be the same as before
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
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
