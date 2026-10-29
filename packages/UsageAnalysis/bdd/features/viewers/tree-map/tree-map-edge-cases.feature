@viewers @edge-cases @realizes:viewers.tree-map
Feature: Tree map edge cases: a split column removed
  Guards for the changes under an open tree map. A split column removed from the table is dropped
  from the split (the viewer's auto-setup keeps only the split names the table has,
  `tree_map_look.dart`). No rows, a renamed split column, a layout naming a missing split column
  and a removed Color or Size column (core#1135) are rows of `viewer-rules-edge-cases.feature`.
  Data: demog-1000 split by RACE, SEX (4 groups, 8 leaves).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tree map viewer with:
      | splitByColumnNames | RACE, SEX |
    Then the "split columns" reading of tree map viewer should be "RACE, SEX"
    And the "leaves" reading of tree map viewer should be 8
    And the "rows shown" reading of tree map viewer should be 1000

  Scenario: A split column removed from the table is dropped from the split
    When user removes "SEX" column
    Then the "split columns" reading of tree map viewer should be "RACE"
    And the "leaves" reading of tree map viewer should be 4
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And no errors should have been logged
