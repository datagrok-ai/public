@viewers @realizes:viewers.pivot-table
Feature: Pivot table regressions from closed issues
  Closed pivot table issues, each claimed where the bug showed, on demog-1000 with the automatic
  cross tab (DIS_POP by SEVERITY, avg(AGE)). The inner grid's colour coding and a column hidden in it
  stay through a filter on the source at Row Source Filtered, and through a selection change at Row
  Source Selected (public#2606). After ADD has published the aggregation as a table, Link Tables
  links the source with it and logs no error (public#3470): row to filter on DIS_POP, so the
  aggregation keeps the one row of the source's current DIS_POP.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    Then the "aggregated rows" reading of pivot table viewer should be 6

  # public#2606
  @realizes:public#2606
  Scenario: At Row Source Filtered the inner grid keeps its colouring and its hidden column through a filter (public#2606)
    When user sets "Row Source" property of pivot table viewer to "Filtered"
    Then the "rows shown" reading of pivot table viewer should be 1000
    When user picks "Grid > Color Coding > Linear" from the context menu of the "grid header None avg(AGE)" area of pivot table viewer
    Then the "color coding of None avg(AGE)" reading of pivot table viewer should be "Linear"
    When user picks "Grid > Hide" from the context menu of the "grid header High avg(AGE)" area of pivot table viewer
    Then the "column visible of High avg(AGE)" reading of pivot table viewer should be "false"
    When user adds a categorical filter on "SEX" keeping "M"
    Then the "rows shown" reading of pivot table viewer should be 447
    And the "color coding of None avg(AGE)" reading of pivot table viewer should be "Linear"
    And the "color of grid cell 2 of None avg(AGE)" reading of pivot table viewer should not be "#ffffff"
    And the "column visible of High avg(AGE)" reading of pivot table viewer should be "false"
    And no errors should have been logged

  # public#2606
  @realizes:public#2606
  Scenario: At Row Source Selected the inner grid keeps its colouring and its hidden column through a selection (public#2606)
    Given user selects rows where "SEX" is "F"
    When user sets "Row Source" property of pivot table viewer to "Selected"
    Then the "rows shown" reading of pivot table viewer should be 553
    When user picks "Grid > Color Coding > Linear" from the context menu of the "grid header None avg(AGE)" area of pivot table viewer
    Then the "color coding of None avg(AGE)" reading of pivot table viewer should be "Linear"
    When user picks "Grid > Hide" from the context menu of the "grid header High avg(AGE)" area of pivot table viewer
    Then the "column visible of High avg(AGE)" reading of pivot table viewer should be "false"
    When user selects rows where "SEX" is "M"
    Then the "rows shown" reading of pivot table viewer should be 447
    And the "color coding of None avg(AGE)" reading of pivot table viewer should be "Linear"
    And the "color of grid cell 2 of None avg(AGE)" reading of pivot table viewer should not be "#ffffff"
    And the "column visible of High avg(AGE)" reading of pivot table viewer should be "false"
    And no errors should have been logged

  # public#3470
  @realizes:public#3470
  Scenario: Link Tables links the source with the table the pivot published (public#3470)
    When user clicks on the "add to workspace" area of pivot table viewer
    Then table "demog-1000 aggregation" should be open
    When user picks "Data > Link Tables..." from the top menu
    And user sets the tables of the Link Tables dialog to "demog-1000" and "demog-1000 aggregation"
    And user sets key columns 1 of the Link Tables dialog to "DIS_POP" and "DIS_POP"
    And user clicks on LINK button in "Link Tables" dialog
    Then UNLINK button in "Link Tables" dialog should be visible
    And 1 rows of table "demog-1000 aggregation" should pass the filter
    When user clicks on CLOSE button in "Link Tables" dialog
    Then no errors should have been logged
    And no error or warning balloon should have been shown
