@viewers @edge-cases @realizes:charts.viewer.tree
Feature: Tree edge cases: renamed and removed columns, empty values, no rows, out-of-range style values
  Guards for the data shapes and the changes under an open Tree. A hierarchy column renamed while
  the tree is open is followed by the hierarchy, also while fewer than 500 rows pass the filter,
  where the viewer skips its category check; the tree keeps drawing and logs no error. A hierarchy
  column removed from the table is dropped from the hierarchy (the viewer's column-removed handler).
  Removing the Size or the Color column throws "No such column" on every later render (the viewer
  drops a removed column from the hierarchy only): two known failures with no ticket yet
  (`@ticket-pending`). An empty value is a branch named `(empty)` only while Include Nulls is on
  (the property's meaning; Charts tree/CLAUDE.md). A filter that removes every row logs no error and
  the branches come back with the rows. Initial Tree Depth declares 0..5 and Font Size at most 30: a
  larger value typed in the context panel logs no error and the tree keeps drawing (whether the
  panel should take such a value is a question in MISSING.md). Data: demog-1000 with the hierarchy
  CONTROL, SEX, RACE.

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

  Scenario: A hierarchy column renamed with every row passing is followed
    When user renames "RACE" column to "ORIGIN"
    Then "hierarchyColumnNames" property of tree viewer should be "CONTROL, SEX, ORIGIN"
    And tree viewer should have a "branch All | false | F | Asian" area
    And no errors should have been logged

  Scenario: A hierarchy column removed from the table is dropped from the hierarchy
    When user removes "RACE" column
    Then "hierarchyColumnNames" property of tree viewer should be "CONTROL, SEX"
    And tree viewer should have a "branch All | false | F" area
    And tree viewer should not have a "branch All | false | F | Asian" area
    And no errors should have been logged

  @known-failure:55 @ticket-pending
  Scenario: Removing the Size column logs no error
    The Tree drops a removed column from its hierarchy only: a removed Size column stays in the
    settings and every later render throws "No such column". No ticket yet. Remove the tag when the
    removal logs no error.
    When user sets "sizeColumnName" property of tree viewer to "HEIGHT"
    And user removes "HEIGHT" column
    Then the table should not have a column "HEIGHT"
    And no errors should have been logged
    And tree viewer should have a "branch All | false | F | Asian" area

  @known-failure:66 @ticket-pending
  Scenario: Removing the Color column logs no error
    The Tree drops a removed column from its hierarchy only: a removed Color column stays in the
    settings and every later render throws "No such column". No ticket yet. Remove the tag when the
    removal logs no error.
    When user sets "colorColumnName" property of tree viewer to "AGE"
    And user removes "AGE" column
    Then the table should not have a column "AGE"
    And no errors should have been logged
    And tree viewer should have a "branch All | false | F | Asian" area

  Scenario: A filter that removes every row logs no error and the branches come back
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And the "rows shown" reading of tree viewer should be 0
    And tree viewer should not have a "branch All | false | F | Asian" area
    When user resets the filter
    Then the "rows shown" reading of tree viewer should be 1000
    And tree viewer should have a "branch All | false | F | Asian" area
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

  Scenario: Values above the declared maximum typed in the context panel log no error
    When user clicks on grid
    And user clicks on settings icon of tree viewer
    Given "Style" category in context panel is expanded
    When user enters "6" into "Initial Tree Depth" property in context panel
    And user enters "31" into "Font Size" property in context panel
    Then tree viewer should have finished updating
    And tree viewer should have a "branch All | false | F | Asian" area
    And no errors should have been logged
