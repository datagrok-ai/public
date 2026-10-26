@viewers @edge-cases @realizes:viewers.tree-map
Feature: Tree map edge cases: no rows, split, color and size columns removed or renamed, a layout with a missing column
  Guards for the data shapes and the changes under an open tree map. A filter that removes every row
  leaves no row on the map and no error, and the leaves come back with the rows. A split column
  removed from the table is dropped from the split (the viewer's auto-setup keeps only the split
  names the table has, `tree_map_look.dart`), and a renamed one is followed (the viewer base renames
  column settings on the table's rename event). A layout that names a split column the table no
  longer has is split by the columns that remain. Removing the Color or the Size column while the
  map is open throws "Column can't be empty." on the removal and on every later refresh (core#1135):
  both scenarios are known failures, and the tags go when core#1135 is fixed. Data: demog-1000 split
  by RACE, SEX (4 groups, 8 leaves).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tree map viewer with:
      | splitByColumnNames | RACE, SEX |
    Then the "split columns" reading of tree map viewer should be "RACE, SEX"
    And the "leaves" reading of tree map viewer should be 8
    And the "rows shown" reading of tree map viewer should be 1000

  Scenario: A filter that removes every row leaves the map empty without an error
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And the "rows shown" reading of tree map viewer should be 0
    And tree map viewer should not have a "leaf Caucasian | M" area
    When user resets the filter
    Then the "rows shown" reading of tree map viewer should be 1000
    And the "leaves" reading of tree map viewer should be 8
    And no errors should have been logged

  Scenario: A split column removed from the table is dropped from the split
    When user removes "SEX" column
    Then the "split columns" reading of tree map viewer should be "RACE"
    And the "leaves" reading of tree map viewer should be 4
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And no errors should have been logged

  Scenario: A split column renamed is followed by the split
    When user renames "SEX" column to "GENDER"
    Then the "split columns" reading of tree map viewer should be "RACE, GENDER"
    And the "leaves" reading of tree map viewer should be 8
    And no errors should have been logged

  Scenario: A layout naming a split column the table no longer has is split by the columns that remain
    When user saves the layout of the current table view
    And user removes "SEX" column
    And user loads the saved layout
    Then the "split columns" reading of tree map viewer should be "RACE"
    And the "leaves" reading of tree map viewer should be 4
    And no errors should have been logged

  @known-failure:60 @realizes:core#1135
  Scenario: Removing the Color column keeps the map working without an error
    core#1135: removing the color column throws "Column can't be empty." on the removal and on every
    later refresh. Remove the tag when core#1135 is fixed.
    When user sets "colorColumnName" property of tree map viewer to "AGE"
    Then the "color column" reading of tree map viewer should be "AGE"
    When user removes "AGE" column
    Then no errors should have been logged
    And the "leaves" reading of tree map viewer should be 8

  @known-failure:70 @realizes:core#1135
  Scenario: Removing the Size column keeps the map working without an error
    core#1135: removing the size column throws "Column can't be empty." on the removal and on every
    later refresh. Remove the tag when core#1135 is fixed.
    When user sets "sizeColumnName" property of tree map viewer to "WEIGHT"
    Then the "size column" reading of tree map viewer should be "WEIGHT"
    When user removes "WEIGHT" column
    Then no errors should have been logged
    And the "leaves" reading of tree map viewer should be 8
