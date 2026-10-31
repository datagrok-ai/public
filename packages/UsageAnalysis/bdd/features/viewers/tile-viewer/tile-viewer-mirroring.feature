@journey @viewers @realizes:viewers.tile-viewer
Feature: Tile viewer mirrors the frame
  A card is a live view of its row, not a copy taken when the lane was built: a cell edited in the
  table reaches the field, a column renamed carries its field, its caption and its value across,
  and a field always shows the grid's display string — the formatted text a cell renders, never the
  raw number behind it — a date format picked for the column in the grid included (public#3855).
  A caption too long for the card shows the whole column name in a tooltip, and a short one shows
  none (public#1593, on a two-row table of its own). This journey mutates the frame, so it owns its
  restore: every scenario puts the cell, the name or the column back (STARTED, given a format by the
  first scenario, is removed and added again by the calculated-column one).
  One journey on demog-1000; row 1 is X0273T21000300003, 26, F, Caucasian, HEIGHT 174.705,
  WEIGHT 74.1, CONTROL false, started 8/2/1990, and row 3 is the first with CONTROL true.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer
    Then tile viewer should be visible
    And the "fields shown" reading of tile viewer should be 10
    And the "tiles" reading of tile viewer should be at least 3

  @realizes:public#3855
  Scenario: A field shows the grid's display string, never the raw cell, a date format picked in the grid included (public#3855)
    Then the "WEIGHT of row 1" reading of tile viewer should be "74.10"
    And the "WEIGHT of row 1" reading of tile viewer should not be "74.1"
    And the "HEIGHT of row 1" reading of tile viewer should be "174.705"
    And the "STARTED of row 1" reading of tile viewer should be "8/2/1990"
    When user scrolls the grid to the "STARTED" column
    And user picks "Format > dd.MM.yyyy (02.08.1990)" from the context menu of the "header STARTED" area of grid
    And user closes the context menu
    Then the "text of cell 1 of STARTED" reading of grid should be "02.08.1990"
    And the "STARTED of row 1" reading of tile viewer should be "02.08.1990"
    And the "CONTROL of row 1" reading of tile viewer should be "false"
    And the "CONTROL of row 3" reading of tile viewer should be "true"
    And no errors should have been logged

  Scenario: A cell edited in the table reaches the card
    Then the "AGE of row 1" reading of tile viewer should be "26"
    When user sets "AGE" column in row 1 to "99"
    Then the value of "AGE" column in row 1 should be "99"
    And the "AGE of row 1" reading of tile viewer should be "99"
    And the "AGE of row 2" reading of tile viewer should be "30"
    When user sets "AGE" column in row 1 to "26"
    Then the "AGE of row 1" reading of tile viewer should be "26"
    And no errors should have been logged

  Scenario: A column rename carries the field, its caption and its value across
    Then tile viewer should have a "field AGE of row 1" area
    And tile viewer should have a "label AGE of row 1" area
    And the "fields" reading of tile viewer should contain "AGE"
    When user renames "AGE" column to "AGE_YRS"
    Then the table should have a column "AGE_YRS"
    And the table should not have a column "AGE"
    And the "fields" reading of tile viewer should contain "AGE_YRS"
    And the "fields" reading of tile viewer should not contain "AGE"
    And tile viewer should have a "field AGE_YRS of row 1" area
    And tile viewer should have a "label AGE_YRS of row 1" area
    And tile viewer should not have a "field AGE of row 1" area
    And the "AGE_YRS of row 1" reading of tile viewer should be "26"
    And the "HEIGHT of row 1" reading of tile viewer should be "174.705"
    And the "fields shown" reading of tile viewer should be 10
    When user renames "AGE_YRS" column to "AGE"
    Then the "fields" reading of tile viewer should contain "AGE"
    And the "AGE of row 1" reading of tile viewer should be "26"
    And tile viewer should have a "field AGE of row 1" area
    And no errors should have been logged

  Scenario: A calculated column that reaches the card shows its formatted text, not the raw number
    When user removes "DEMOG" column
    And user removes "STARTED" column
    Then the table should have 9 columns
    When user adds a calculated column "H_THIRD" with formula "${HEIGHT} / 3"
    Then the table should have 10 columns
    And the "fields shown" reading of tile viewer should be 10
    And the "fields" reading of tile viewer should contain "H_THIRD"
    And tile viewer should have a "field H_THIRD of row 1" area
    And the "H_THIRD of row 1" reading of tile viewer should be "58.24"
    And the "H_THIRD of row 1" reading of tile viewer should not be "58.235"
    And no errors should have been logged
    When user removes "H_THIRD" column
    And user adds a calculated column "DEMOG" with formula "${AGE}"
    And user adds a calculated column "STARTED" with formula "${AGE}"
    Then the table should have 11 columns

  @realizes:public#1593
  Scenario: A column name the card cuts off shows whole in a tooltip over its caption (public#1593)
    Given user opens a table "long-names" with:
      | id | measured_concentration_of_compound_in_plasma |
      | 1  | 5                                            |
      | 2  | 6                                            |
    And user adds a tile viewer
    When user hovers over the "label measured_concentration_of_compound_in_plasma of row 1" area of tile viewer
    Then tooltip should be visible
    And tooltip should contain text "measured_concentration_of_compound_in_plasma"
    When user hovers over the "label id of row 1" area of tile viewer
    Then tooltip should be hidden
    When user moves the pointer away from tile viewer
    Then no errors should have been logged
