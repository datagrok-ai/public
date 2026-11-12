@journey @viewers @realizes:viewers.grid
Feature: An int column shows the named number formats
  The named formats of a float column work on an int column too, instead of being read as a number
  pattern, which printed the name in front of the value (core#249). The first AGE of demog-1000 is 26.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  @realizes:core#249
  Scenario: An int column takes the float format names
    When user sets the format of "AGE" column in the grid to "scientific"
    Then the "text of cell 1 of AGE" reading of grid should be "2.60E1"
    When user sets the format of "AGE" column in the grid to "two digits after comma"
    Then the "text of cell 1 of AGE" reading of grid should be "26.00"
    And no errors should have been logged
