@journey @viewers @realizes:viewers.stats-viewer
Feature: Which rows the statistics are computed over, and which columns get a row
  The statistics are computed over `combinedFilter`, so **Row Source** and the table's own filter
  both move them. The old spec set Row Source, applied a filter, and asserted a five-hundred-pixel
  canvas delta — which says something was repainted and nothing about which numbers changed. Here
  the claim is the numbers: 1000 rows average 45.68, the 447 male rows average 44.55 and the 553
  female rows average 46.59, and `rows shown` and `values of AGE` say how many rows were counted.
  The **Columns** scenario replaces an old assertion that read the property grid's column-picker
  caption for the text "/ 11" — a claim about the picker widget, not about this viewer. `columns`
  is what the table has and `columns shown` is how many got a row, and a column that got none has
  no `row <COLUMN>` area.
  The Background opens the filter panel before it sizes the viewer, and that is not decoration: a
  cell is reported only for a column the grid scrolled into view, so a filter panel opening
  half-way through the journey narrows the viewer and takes `avg` off screen — and with it the
  `avg of AGE` reading. Fix the geometry first and the readings stay comparable across the whole
  feature.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a statistics viewer
    And user resizes statistics viewer to 900 by 400
    Then 1000 rows should pass the filter
    And "rowSource" property of statistics viewer should be "Filtered"
    And the "rows shown" reading of statistics viewer should be 1000
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    And statistics viewer should report no error

  Scenario: With Row Source Filtered the table's filter moves every number
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of statistics viewer should be 447
    And the "values of AGE" reading of statistics viewer should be "447"
    And the "avg of AGE" reading of statistics viewer should be "44.55"
    And the "unique of SEX" reading of statistics viewer should be "1"
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows shown" reading of statistics viewer should be 1000
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    And the "unique of SEX" reading of statistics viewer should be "2"
    And no errors should have been logged

  Scenario: With Row Source Selected the selection is what is counted
    When user sets "rowSource" property of statistics viewer to "Selected"
    And user selects rows where "SEX" is "F"
    Then 553 rows should be selected
    And the "rows shown" reading of statistics viewer should be 553
    And the "values of AGE" reading of statistics viewer should be "553"
    And the "avg of AGE" reading of statistics viewer should be "46.59"
    And the "max of AGE" reading of statistics viewer should be "80.00"
    When user selects rows where "SEX" is "M"
    Then the "rows shown" reading of statistics viewer should be 447
    And the "avg of AGE" reading of statistics viewer should be "44.55"
    And the "max of AGE" reading of statistics viewer should be "89.00"
    When user selects no rows
    And user sets "rowSource" property of statistics viewer to "All"
    Then the "rows shown" reading of statistics viewer should be 1000
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    When user sets "rowSource" property of statistics viewer to "Filtered"
    Then no errors should have been logged

  Scenario: With Row Source All the table's filter leaves the statistics where they are
    When user sets "rowSource" property of statistics viewer to "All"
    And user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of statistics viewer should be 1000
    And the "values of AGE" reading of statistics viewer should be "1000"
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    When user sets "rowSource" property of statistics viewer to "Filtered"
    Then the "rows shown" reading of statistics viewer should be 447
    And the "avg of AGE" reading of statistics viewer should be "44.55"
    When user resets the filter
    Then the "rows shown" reading of statistics viewer should be 1000
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    And no errors should have been logged

  Scenario: Columns decides which of the table's columns gets a row
    Then the "columns" reading of statistics viewer should be 11
    And the "columns shown" reading of statistics viewer should be 11
    And statistics viewer should have a "row HEIGHT" area
    When user sets "columnNames" property of statistics viewer to "USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    Then the "columns shown" reading of statistics viewer should be 10
    And the "columns" reading of statistics viewer should be 11
    And statistics viewer should not have a "row HEIGHT" area
    And statistics viewer should have a "row WEIGHT" area
    When user sets "columnNames" property of statistics viewer to "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    Then the "columns shown" reading of statistics viewer should be 11
    And statistics viewer should have a "row HEIGHT" area
    And the "values of HEIGHT" reading of statistics viewer should be "872"
    And no errors should have been logged
