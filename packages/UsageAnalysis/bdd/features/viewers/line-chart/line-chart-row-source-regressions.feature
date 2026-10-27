@viewers @realizes:viewers.line-chart
Feature: Line chart regressions — row source, row groups, categories and the log axis
  Closed line chart issues, each claimed where the bug showed: a layout with a categorical X and Row
  Source = Selected, applied with nothing selected, loads without an error (public#2350); a split
  chart switched to Row Source = Selected logs nothing and then draws exactly the selected rows
  (public#1650, public#1628); a row group hovered in a bar chart repaints the line chart only while
  Show Mouse Over Row Group is on (public#3709); a molecule column on X draws a category per
  structure (public#1582); and a logarithmic Y axis fits the data instead of leaving an empty band
  above it (public#2485). Each scenario adds its own line chart. The split legend under Row Source =
  Selected is `line-chart-regressions.feature`.
  Not translated: the line click under a split (public#2608, public#2360), on-plot split and Y
  column removal (public#2951 and its duplicates), the log-axis warning (public#3701), label
  alignment after a filter (public#1456), a custom category order (public#1551) and the other
  candidates listed in MISSING.md — each needs a step or a core signal that does not exist yet.

  Background:
    Given user is logged in

  # public#2350
  @realizes:public#2350
  Scenario: A layout with a categorical X and Row Source Selected loads with nothing selected (public#2350)
    Given user opens demog-1000 dataset
    When user adds a line chart viewer with:
      | xColumnName  | RACE     |
      | yColumnNames | AGE      |
      | rowSource    | Selected |
    And user selects rows where "RACE" is "Asian"
    Then the "rows shown" reading of line chart viewer should be 15
    When user saves the layout of the current table view
    And user clears the row selection
    Then the "rows shown" reading of line chart viewer should be 0
    When user loads the saved layout
    Then line chart viewer should report no error
    And "rowSource" property of line chart viewer should be "Selected"
    And the "rows shown" reading of line chart viewer should be 0
    And no errors should have been logged
    When user selects rows where "RACE" is "Black"
    Then the "rows shown" reading of line chart viewer should be 27
    And the "x categories" reading of line chart viewer should be 1
    And no errors should have been logged

  # public#1650, public#1628: Pyrrolidines holds 21 rows at 20 distinct masses (the chart aggregates
  # the pair), Aminopiperidines 4 at 4
  @realizes:public#1650 @realizes:public#1628
  Scenario: A split chart under Row Source Selected logs nothing and draws exactly the selected rows (public#1650)
    Given user opens spgi dataset
    When user adds a line chart viewer with:
      | xColumnName      | Average Mass     |
      | yColumnNames     | Chemical Space X |
      | splitColumnNames | Series           |
    Then the "lines" reading of line chart viewer should be 5
    When user sets "rowSource" property of line chart viewer to "Selected"
    Then the "rows shown" reading of line chart viewer should be 0
    And line chart viewer should report no error
    And no errors should have been logged
    When user selects rows where "Series" is "Pyrrolidines"
    Then line chart viewer should show the selected rows of its table
    And the "rows shown" reading of line chart viewer should be 21
    And the "lines" reading of line chart viewer should be 1
    And the "markers drawn" reading of line chart viewer should be 20
    When user selects rows where "Series" is one of "Pyrrolidines, Aminopiperidines"
    Then the "rows shown" reading of line chart viewer should be 25
    And the "lines" reading of line chart viewer should be 2
    And the "markers drawn" reading of line chart viewer should be 24
    And no errors should have been logged

  # public#3709
  @realizes:public#3709
  Scenario: A row group hovered in a bar chart repaints the line chart only while Show Mouse Over Row Group is on (public#3709)
    Given user opens demog-1000 dataset
    When user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
    And user adds a bar chart viewer with:
      | Split | RACE |
    And user takes a snapshot of line chart viewer
    And user hovers over the "bar Asian" area of bar chart viewer
    Then line chart viewer should have repainted
    When user moves the pointer away from grid
    And user sets "showMouseOverRowGroup" property of line chart viewer to "false"
    And user takes a snapshot of line chart viewer
    And user hovers over the "bar Asian" area of bar chart viewer
    Then bar chart viewer should have repainted
    And line chart viewer should not have repainted
    When user moves the pointer away from grid
    And user sets "showMouseOverRowGroup" property of line chart viewer to "true"
    Then no errors should have been logged

  # public#1582
  @realizes:public#1582
  Scenario: A molecule column on X draws a category per structure (public#1582)
    Given user opens spgi dataset
    When user adds a line chart viewer with:
      | xColumnName  | Core             |
      | yColumnNames | Chemical Space X |
    Then the "x categories" reading of line chart viewer should be 6
    And line chart viewer should report no error
    And line chart viewer should be painted
    And no errors should have been logged

  # public#2485: the column's values run from 0.04 to 31.5 (two zeros, which a log axis drops)
  @realizes:public#2485
  Scenario: A logarithmic Y axis fits the data (public#2485)
    Given user opens spgi dataset
    When user adds a line chart viewer with:
      | xColumnName  | Num Heavy Atoms     |
      | yColumnNames | Whole blood assay 1 |
    And user sets "yAxisType" property of line chart viewer to "logarithmic"
    Then the "y axis max of \"Whole blood assay 1\"" reading of line chart viewer should be between 31.5 and 37.8
    And the "y axis min of \"Whole blood assay 1\"" reading of line chart viewer should be between 0.04 and 0.14
    And line chart viewer should report no error
    And no errors should have been logged
