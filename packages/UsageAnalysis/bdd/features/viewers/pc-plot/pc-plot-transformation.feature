@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot transformations
  A Transformation aggregates the table first and the plot draws the aggregated frame: a pivot
  replaces the axes with the pivot's categories, an aggregation with the aggregated columns, and
  clearing it puts the original axes back. A transformation that does not parse raises a balloon and
  leaves the axes alone, and closing the filter panel while one is applied must not break the view
  (GROK-18091).
  spgi-100: 100 rows whose ten numerical columns the plot auto-picks. `Chemist 521` has 13 distinct
  values and `Series` five (Triazoles, Diazabicyclooctane, Pyrrolidines, the blank one and
  Aminopiperidines), so keying on the chemist and pivoting on the series gives 13 rows and those
  five axes; keying on the series and averaging two numbers gives 5 rows and two axes.
  Known open bug: the last scenario is GROK-17306 — with a transformation applied, the filter
  panel's Reset filters drops the row selection instead of leaving it alone. It is translated as
  written and is expected to fail; it runs last so nothing after it inherits its state.
  spgi-100 carries a molecule column, so the filter panel builds a substructure filter for it: with
  the package autostarts still pending the platform logs `Cannot execute "Chem:substructureFilter"
  synchronously, either the package is not loaded or the function is not there`. That is the panel's
  own precondition and nothing to do with the PC plot, so the scenario waits for the autostarts the
  platform already announces instead of tolerating the error — a filter panel opened on a molecule
  table before its package is up is a finding of its own.
  Defect found while writing this: `_invalidTransformation` is only cleared when a later
  transformation parses, so after an invalid one is cleared to "" the plot keeps reporting
  `error` = "Invalid transformation" — the scenario below restores a valid transformation to clear
  it, and states that as the product's behaviour rather than papering over it.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a pc plot viewer
    Then the "transformed" reading of pc plot viewer should be "false"
    And the "axes" reading of pc plot viewer should be 10
    And pc plot viewer should show 100 rows
    And the "error" reading of pc plot viewer should be ""
    And user remembers the "axis order" reading of pc plot viewer

  Scenario: A pivot replaces the axes with the pivot column's categories
    When user sets "Transformation" property of pc plot viewer to '[{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]'
    Then the "transformed" reading of pc plot viewer should be "true"
    And the "axes" reading of pc plot viewer should be 5
    And the axes of pc plot viewer should be "Triazoles, Diazabicyclooctane, Pyrrolidines, , Aminopiperidines"
    And pc plot viewer should have an "axis \"Triazoles\"" area
    And the "rows" reading of pc plot viewer should be 13
    And pc plot viewer should show 13 rows
    And the "axis order" reading of pc plot viewer should not be as remembered
    And pc plot viewer should be painted
    And no errors should have been logged

  Scenario: Clearing the transformation puts the original axes back
    When user sets "Transformation" property of pc plot viewer to ""
    Then the "transformed" reading of pc plot viewer should be "false"
    And the "axes" reading of pc plot viewer should be 10
    And the "axis order" reading of pc plot viewer should be as remembered
    And pc plot viewer should not have an "axis \"Triazoles\"" area
    And pc plot viewer should show 100 rows
    And no errors should have been logged

  Scenario: An aggregation draws the aggregated columns
    When user sets "Transformation" property of pc plot viewer to '[{"#type":"GroupAggregation","aggType":"key","colName":"Series"},{"#type":"GroupAggregation","aggType":"avg","colName":"Average Mass"},{"#type":"GroupAggregation","aggType":"avg","colName":"TPSA"}]'
    Then the "transformed" reading of pc plot viewer should be "true"
    And the axes of pc plot viewer should be "avg(Average Mass), avg(TPSA)"
    And the "rows" reading of pc plot viewer should be 5
    And pc plot viewer should show 5 rows
    And pc plot viewer should be painted
    When user sets "Transformation" property of pc plot viewer to ""
    Then the "transformed" reading of pc plot viewer should be "false"
    And the "axis order" reading of pc plot viewer should be as remembered
    And no errors should have been logged

  Scenario: A transformation that does not parse raises a balloon and leaves the axes alone
    When user sets "Transformation" property of pc plot viewer to "not a transformation"
    Then an error balloon containing "Invalid transformation" should have been shown
    And the "error" reading of pc plot viewer should be "Invalid transformation"
    And the "transformed" reading of pc plot viewer should be "false"
    And the "axis order" reading of pc plot viewer should be as remembered
    And pc plot viewer should show 100 rows
    And pc plot viewer should be painted
    When user sets "Transformation" property of pc plot viewer to '[{"#type":"GroupAggregation","aggType":"key","colName":"Series"},{"#type":"GroupAggregation","aggType":"avg","colName":"Average Mass"},{"#type":"GroupAggregation","aggType":"avg","colName":"TPSA"}]'
    Then the "error" reading of pc plot viewer should be ""
    And the "transformed" reading of pc plot viewer should be "true"
    When user sets "Transformation" property of pc plot viewer to ""
    Then the "axis order" reading of pc plot viewer should be as remembered
    And no errors should have been logged

  Scenario: Closing the filter panel with a transformation applied breaks nothing (GROK-18091)
    Given the package autostarts have completed
    When user opens the filter panel
    Then filter panel should be visible
    When user sets "Transformation" property of pc plot viewer to '[{"#type":"GroupAggregation","aggType":"key","colName":"Series"},{"#type":"GroupAggregation","aggType":"avg","colName":"Average Mass"},{"#type":"GroupAggregation","aggType":"avg","colName":"TPSA"}]'
    Then the "transformed" reading of pc plot viewer should be "true"
    And pc plot viewer should show 5 rows
    When user clicks on close icon of filters viewer
    Then filter panel should be absent
    And pc plot viewer should be visible
    And pc plot viewer should be painted
    And the "transformed" reading of pc plot viewer should be "true"
    And the "rows" reading of pc plot viewer should be 5
    And the table should have 100 rows
    And no error or warning balloon should have been shown
    When user sets "Transformation" property of pc plot viewer to ""
    Then the "axis order" reading of pc plot viewer should be as remembered
    And no errors should have been logged

  @known-failure
  Scenario: With a transformation, Reset filters restores the rows but drops the selection (GROK-17306)
    Given user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then pc plot viewer should show 1000 rows
    When user sets "Transformation" property of pc plot viewer to '[{"#type":"GroupAggregation","aggType":"key","colName":"SEX"},{"#type":"GroupAggregation","aggType":"pivot","colName":"DIS_POP"},{"#type":"GroupAggregation","aggType":"avg","colName":"WEIGHT"}]'
    Then the "transformed" reading of pc plot viewer should be "true"
    When user opens the filter panel
    And user adds a range filter on "AGE" from 30 to 50
    Then 494 rows should pass the filter
    When user selects the first 10 rows
    Then 10 rows should be selected
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And 10 rows should be selected
