@viewers @realizes:viewers.form
Feature: Form viewer edge cases
  Guards for the form viewer at the ends of its table and for the toolbar's Load form from file, on
  demog-1000 with a default form of its 11 columns. The previous-row arrow on the first row and the
  next-row arrow on the last row leave the current row where it is (`form_core.dart` moves only to a
  row inside the table). A file that is not a form, a CSV file or a JSON file with no form in it (a
  saved table layout), is refused with an error balloon "<file> is not a valid form file" and the form
  keeps its fields (core#1138). The files are fixtures the project already keeps
  (`fixtures/cars-small.csv`, `fixtures/nx/spgi-legend-position-old.layout`). In an editable form,
  text typed into a numeric field leaves the cell as it was (core#1192), and a column removed from
  the table takes its field with it, on a default form and on a form built from chosen columns
  (today the second keeps an empty field).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a form viewer
    Then the "fields shown" reading of form viewer should be 11

  @edge-cases
  Scenario: The row arrows stop at the first and the last row
    When user clicks on the "next row" area of form viewer
    Then row 2 should be current
    When user clicks on the "previous row" area of form viewer
    Then row 1 should be current
    When user clicks on the "previous row" area of form viewer
    Then row 1 should be current
    And the "row" reading of form viewer should be 1
    When user makes the last row current
    Then the "row" reading of form viewer should be 1000
    When user clicks on the "next row" area of form viewer
    Then row 1000 should be current
    And the "row" reading of form viewer should be 1000
    And no errors should have been logged

  @edge-cases @realizes:core#1138
  Scenario Outline: A file that is not a form is refused with a message and the form keeps its fields
    When user uploads "<file>" through folder open icon in form viewer
    Then an error balloon containing "<name> is not a valid form file" should have been shown
    And the "fields shown" reading of form viewer should be 11
    And no errors should have been logged

    Examples:
      | file                                        | name                             |
      | fixtures/cars-small.csv                     | cars-small.csv                   |
      | fixtures/nx/spgi-legend-position-old.layout | spgi-legend-position-old.layout  |

  @edge-cases @known-failure:54 @realizes:core#1192
  Scenario: Text typed into a numeric field of an editable form leaves the cell as it was
    When user clicks on the "edit" area of form viewer
    And user enters "31" into the "field AGE" area of form viewer
    Then the value of "AGE" column in row 1 should be "31"
    When user enters "abc" into the "field AGE" area of form viewer
    Then the value of "AGE" column in row 1 should be "31"
    And no errors should have been logged

  @edge-cases
  Scenario: A column removed from the table takes its field off the form
    When user removes "AGE" column
    Then the "fields shown" reading of form viewer should be 10
    And form viewer should not have a "field AGE" area
    And form viewer should report no error
    And no errors should have been logged

  @edge-cases @known-failure:70 @ticket-pending
  Scenario: A column removed from the table takes its field off a form built from chosen columns
    Given user sets "columnNames" property of form viewer to "AGE, HEIGHT, WEIGHT"
    Then the "fields shown" reading of form viewer should be 3
    When user removes "WEIGHT" column
    Then the "fields shown" reading of form viewer should be 2
    And form viewer should not have a "field WEIGHT" area
    And no errors should have been logged
