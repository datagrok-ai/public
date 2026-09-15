@journey @viewers @realizes:viewers.tree-map
Feature: Tree map splitting, nesting and the on-viewer selectors
  Which column the map splits by, how many rectangles that makes, and what each one holds.
  The spec this replaces never named a rectangle: `hoverLeaf()` moved the pointer to five
  hard-coded fractions of the viewer box until a tooltip matched `/(\d+) rows/`, threw when none
  did, and then compared the number it had parsed against a `groupBy` re-implemented inside the
  test. Every node of the last layout is now a hit area addressed by its path — `leaf Caucasian`,
  `group Caucasian`, `leaf Caucasian | M` — and `rows of <path>` is the viewer's own count, so
  there is nothing left to hunt for and nothing left to re-derive.
  The auto split is named rather than asserted non-empty: `TreeMapLook.auto` takes the first
  categorical column with 5..19 categories, which on demog-1000 is DIS_POP (6), not RACE.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tree map viewer
    Then the "rows shown" reading of tree map viewer should be 1000
    And the "error" reading of tree map viewer should be ""
    And tree map viewer should be painted

  Scenario: The map opens split by the first column with 5 to 19 categories
    Then the "split columns" reading of tree map viewer should be "DIS_POP"
    And the "levels" reading of tree map viewer should be 1
    And the "selectors" reading of tree map viewer should be 2
    And the "leaves" reading of tree map viewer should be 6
    And the "groups" reading of tree map viewer should be 0
    And the "rows of RA" reading of tree map viewer should be 434
    And the "rows of Psoriasis" reading of tree map viewer should be 204
    And the "rows of PsA" reading of tree map viewer should be 38
    And tree map viewer should have a "leaf RA" area
    And tree map viewer should have a "leaf body RA" area
    And the "leaf RA" area of tree map viewer should be painted
    And no errors should have been logged

  Scenario: Splitting by RACE gives one rectangle per race holding exactly its rows
    When user sets "splitByColumnNames" property of tree map viewer to "RACE"
    Then the "split columns" reading of tree map viewer should be "RACE"
    And the "leaves" reading of tree map viewer should be 4
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And the "rows of Other" reading of tree map viewer should be 62
    And the "rows of Black" reading of tree map viewer should be 27
    And the "rows of Asian" reading of tree map viewer should be 15
    And the "area of Caucasian" reading of tree map viewer should be 896
    And the "rows shown" reading of tree map viewer should be 1000
    And tree map viewer should have repainted
    And tree map viewer should have a "leaf Caucasian" area
    And tree map viewer should not have a "leaf RA" area
    When user sets "splitByColumnNames" property of tree map viewer to "DIS_POP"
    Then the "leaves" reading of tree map viewer should be 6
    And no errors should have been logged

  Scenario: A second level nests each race into its sexes and the group sums them
    When user sets "splitByColumnNames" property of tree map viewer to "RACE, SEX"
    Then the "levels" reading of tree map viewer should be 2
    And the "selectors" reading of tree map viewer should be 3
    And the "leaves" reading of tree map viewer should be 8
    And the "groups" reading of tree map viewer should be 4
    And the "rows of Caucasian | M" reading of tree map viewer should be 416
    And the "rows of Caucasian | F" reading of tree map viewer should be 480
    And the "rows of Caucasian" reading of tree map viewer should be 896
    And tree map viewer should have a "group Caucasian" area
    And tree map viewer should have a "leaf Caucasian | M" area
    And tree map viewer should not have a "leaf Caucasian" area
    And tree map viewer should have repainted
    When user sets "splitByColumnNames" property of tree map viewer to "DIS_POP"
    Then the "levels" reading of tree map viewer should be 1
    And the "groups" reading of tree map viewer should be 0
    And the "selectors" reading of tree map viewer should be 2
    And no errors should have been logged

  Scenario: Show Column Selection Panel takes the selectors off the map
    Then the "selection panel shown" reading of tree map viewer should be "true"
    And tree map viewer should have a "split selector 1" area
    And tree map viewer should have a "split selector 2" area
    And tree map viewer should have a "color selector" area
    When user sets "showColumnSelectionPanel" property of tree map viewer to "false"
    Then the "selection panel shown" reading of tree map viewer should be "false"
    And tree map viewer should not have a "split selector 1" area
    And tree map viewer should not have a "split selector 2" area
    And tree map viewer should not have a "color selector" area
    And the "levels" reading of tree map viewer should be 1
    When user sets "showColumnSelectionPanel" property of tree map viewer to "true"
    Then the "selection panel shown" reading of tree map viewer should be "true"
    And tree map viewer should have a "split selector 1" area
    And no errors should have been logged

  Scenario: The title bar closes the map
    When user clicks on close icon of tree map viewer
    Then tree map viewer should be absent
    And the open tableview should have 0 tree map viewers
    And no errors should have been logged
