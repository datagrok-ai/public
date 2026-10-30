@viewers @edge-cases @realizes:charts.viewer.tree
Feature: Tree edge cases: a column renamed under a filter, a hierarchy column removed, empty values, On Click switched from Filter to Select
  Guards for the data shapes and the changes under an open Tree. A hierarchy column renamed while
  fewer than 500 rows pass the filter, where the viewer skips its category check, is followed by the
  hierarchy; the tree keeps drawing and logs no error. A hierarchy column removed from the table is
  dropped from the hierarchy (the viewer's column-removed handler). An empty value is a branch named
  `(empty)` only while Include Nulls is on (the property's meaning; Charts tree/CLAUDE.md).
  Switching On Click from Filter to Select gives back the rows a clicked branch filtered out
  (public#4425). No rows, a rename with every row passing, a removed
  Size or Color column (public#4423), a layout naming a missing column, values typed out of bounds
  and the filter left by a closed tree are rows of `viewer-rules-edge-cases.feature`. Data:
  demog-1000 with the hierarchy CONTROL, SEX, RACE.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    And user adds a tree viewer with:
      | hierarchyColumnNames | CONTROL, SEX, RACE |
    Then tree viewer should have a "branch All | false | F | Asian" area
    And the "rows shown" reading of tree viewer should be 1000

  Scenario: A hierarchy column renamed under a filter of fewer than 500 rows is followed
    When user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    When user renames "RACE" column to "ORIGIN"
    Then "hierarchyColumnNames" property of tree viewer should be "CONTROL, SEX, ORIGIN"
    And tree viewer should have a "branch All | false | M | Asian" area
    And the "rows shown" reading of tree viewer should be 447
    And no errors should have been logged

  Scenario: A hierarchy column removed from the table is dropped from the hierarchy
    When user removes "RACE" column
    Then "hierarchyColumnNames" property of tree viewer should be "CONTROL, SEX"
    And tree viewer should have a "branch All | false | F" area
    And tree viewer should not have a "branch All | false | F | Asian" area
    And no errors should have been logged

  Scenario: An empty value is the (empty) branch only while Include Nulls is on
    Given user opens a table "with-empty" with:
      | site | arm |
      | A    | x   |
      | A    |     |
      | B    | y   |
    And user adds a tree viewer with:
      | hierarchyColumnNames | site, arm |
      | includeNulls         | true      |
    Then tree viewer should have a "branch All | A | (empty)" area
    When user sets "includeNulls" property of tree viewer to "false"
    Then tree viewer should not have a "branch All | A | (empty)" area
    And tree viewer should have a "branch All | A | x" area
    When user sets "includeNulls" property of tree viewer to "true"
    Then tree viewer should have a "branch All | A | (empty)" area
    And no errors should have been logged

  @realizes:public#4425
  Scenario: Switching On Click from Filter to Select gives back the rows its clicked branch filtered out
    When user sets "onClick" property of tree viewer to "Filter"
    And user clicks on the "branch All | false | F | Asian" area of tree viewer
    Then 7 rows should pass the filter
    When user sets "onClick" property of tree viewer to "Select"
    Then all rows should pass the filter
    And no errors should have been logged
