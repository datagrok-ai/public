@viewers @realizes:viewers.form
Feature: Form viewer regressions from closed issues
  Closed form viewer issues, each claimed where the bug showed. An empty cell of an integer column
  shows an empty field, not the integer null -2147483648 (public#3561), on a three-row inline table.
  On demog-1000 sorted in the grid, a field's colour coding is the colour of the row the form shows,
  not of the row that sits at that place in the grid (public#940), and the form's next-row arrow
  takes the next row of the grid's order, not of the table's (public#2644): sorted by HEIGHT
  descending, row 1 (174.705) is followed by row 808 (174.652), not by row 2.

  Background:
    Given user is logged in

  # public#3561
  @realizes:public#3561
  Scenario: An empty integer cell shows an empty field on the form viewer (public#3561)
    Given user opens a table "ints" with:
      | id | n |
      | 1  | 5 |
      | 2  |   |
      | 3  | 7 |
    Then "n" column should have type "int"
    When user adds a form viewer
    And user makes row 1 current
    Then the "n" reading of form viewer should be "5"
    When user makes row 2 current
    Then the "row" reading of form viewer should be 2
    And the "n" reading of form viewer should be ""
    And no errors should have been logged

  # public#940
  @realizes:public#940
  Scenario: A field keeps its own row's colour after the grid is sorted (public#940)
    Given user opens demog-1000 dataset
    And user adds a form viewer
    When user colors "AGE" column linearly from "#FF0000" to "#0000FF"
    And user makes row 5 current
    Then the "background of AGE" reading of form viewer should not be ""
    When user remembers the "background of AGE" reading of form viewer
    And user makes row 7 current
    Then the "background of AGE" reading of form viewer should not be as remembered
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    When user makes row 5 current
    Then the "row" reading of form viewer should be 5
    And the "background of AGE" reading of form viewer should be as remembered
    And no errors should have been logged

  # public#2644
  @realizes:public#2644
  Scenario: After the grid is sorted the next-row arrow takes the next row of the grid (public#2644)
    Given user opens demog-1000 dataset
    And user adds a form viewer
    When user double-clicks on the "header HEIGHT" area of grid
    Then the "sort column" reading of grid should be "HEIGHT"
    And the "sort direction" reading of grid should be "descending"
    When user makes row 1 current
    Then the "HEIGHT" reading of form viewer should be "174.705"
    When user clicks on the "next row" area of form viewer
    Then row 808 should be current
    And the "row" reading of form viewer should be 808
    And the "HEIGHT" reading of form viewer should be "174.652"
    And no errors should have been logged
