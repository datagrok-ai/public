@viewers @edge-cases @realizes:viewers.tile-viewer
Feature: Tile viewer edge cases: nothing selected under Row Source Selected, one row
  Guards for the data shapes under an open tile viewer, which the platform's generic
  data-manipulation test does not reach. Row Source Selected with nothing selected leaves no card
  and no error, and the cards come with the selection. A one-row table shows its one card. No rows,
  the lanes column removed or renamed, and an empty lane value are rows of
  `viewer-rules-edge-cases.feature`. Data: demog-1000 with lanes on RACE (Asian, Black, Caucasian,
  Other).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer with:
      | Lanes Column Name | RACE |
    Then tile viewer should show 1000 rows
    And the "lane names" reading of tile viewer should be "Asian, Black, Caucasian, Other"

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

  Scenario: Escape in the form designer, a view and not a dialog, leaves it open with the field deleted
    When user picks "Edit Form..." from the viewer menu of tile viewer
    Then form designer should be visible
    When user deletes the "AGE" value field in the form designer
    Then the form designer should show value fields "CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    When user presses Escape
    Then form designer should be visible
    And the form designer should show value fields "CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    And no errors should have been logged

  @known-failure:62 @realizes:core#1217
  Scenario: A card dropped into the (empty) lane gets an empty value
    Today the lane of the empty value has no name (the empty-category row of
    `viewer-rules-edge-cases.feature`), so there is no "(empty)" lane to drop the card into.
    Given user opens a table "with-empty" with:
      | site | arm |
      | A    | x   |
      | A    |     |
      | B    | y   |
    And user adds a tile viewer with:
      | Lanes Column Name | arm |
    Then the "lane of row 1" reading of tile viewer should be "x"
    When user moves the pointer away from tile viewer
    And user drags the card of row 1 of tile viewer into lane "(empty)"
    Then the value of "arm" column in row 1 should be ""
    And the "lane of row 1" reading of tile viewer should be "(empty)"
    And no errors should have been logged

  Scenario: A card dragged between the lanes of a calculated column writes the lane's value
    The edit is the user's: AGE_GROUP takes the value of the lane and AGE stays as it was.
    When user adds a calculated column "AGE_GROUP" with formula "if(${AGE} > 40, \"old\", \"young\")"
    And user sets "lanesColumnName" property of tile viewer to "AGE_GROUP"
    Then the "lane of row 1" reading of tile viewer should be "young"
    When user moves the pointer away from tile viewer
    And user drags the card of row 1 of tile viewer into lane "old"
    Then the "lane of row 1" reading of tile viewer should be "old"
    And the value of "AGE_GROUP" column in row 1 should be "old"
    And the value of "AGE" column in row 1 should be "26"
    And no errors should have been logged
