@viewers @edge-cases @realizes:viewers.tile-viewer
Feature: Tile viewer edge cases: no rows, one row, the lanes column removed or renamed
  Guards for the data shapes and the changes under an open tile viewer, which the platform's generic
  data-manipulation test does not reach. A filter that removes every row, and Row Source Selected
  with nothing selected, leave no card and no error; the cards come back with the rows. A one-row
  table shows its one card. The lanes column renamed is followed by the viewer's settings (the
  viewer base renames column settings on the table's rename event), so the lanes stay as they were.
  The lanes column removed while the lanes are shown logs no error and the viewer keeps every row;
  which lanes it should show then is not defined: see MISSING.md. Data: demog-1000 with lanes on
  RACE (Asian, Black, Caucasian, Other).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer with:
      | Lanes Column Name | RACE |
    Then tile viewer should show 1000 rows
    And the "lane names" reading of tile viewer should be "Asian, Black, Caucasian, Other"

  Scenario: A filter that removes every row leaves no card and no error
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And tile viewer should show 0 rows
    And the "tiles" reading of tile viewer should be 0
    And tile viewer should not have a "tile of row 1" area
    When user resets the filter
    Then tile viewer should show 1000 rows
    And tile viewer should have a "tile of row 1" area
    And no errors should have been logged

  Scenario: Row Source Selected with nothing selected shows no card and no error
    Given user selects no rows
    When user sets "Row Source" property of tile viewer to "Selected"
    Then tile viewer should show 0 rows
    And the "tiles" reading of tile viewer should be 0
    When user selects rows where "RACE" is "Asian"
    Then tile viewer should show 15 rows
    And the "lane of row 10" reading of tile viewer should be "Asian"
    And no errors should have been logged

  Scenario: A one-row table shows its one card
    Given user opens a table "one-row" with:
      | name | group | x  |
      | a    | g1    | 42 |
    And user adds a tile viewer with:
      | Lanes Column Name | group |
    Then tile viewer should show 1 rows
    And the "tiles" reading of tile viewer should be 1
    And the "lane of row 1" reading of tile viewer should be "g1"
    And no errors should have been logged

  Scenario: The lanes column renamed is followed and the lanes stay
    When user renames "RACE" column to "ORIGIN"
    Then "Lanes Column Name" property of tile viewer should be "ORIGIN"
    And the "lane names" reading of tile viewer should be "Asian, Black, Caucasian, Other"
    And the "lane of row 10" reading of tile viewer should be "Asian"
    And no errors should have been logged

  Scenario: The lanes column removed while the lanes are shown logs no error
    When user removes "RACE" column
    Then the table should not have a column "RACE"
    And tile viewer should show 1000 rows
    And no errors should have been logged
