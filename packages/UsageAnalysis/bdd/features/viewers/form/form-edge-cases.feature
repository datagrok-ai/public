@viewers @realizes:viewers.form
Feature: Form viewer edge cases
  Guards for the form viewer at the ends of its table and for the toolbar's Load form from file, on
  demog-1000 with a default form of its 11 columns. The previous-row arrow on the first row and the
  next-row arrow on the last row leave the current row where it is (`form_core.dart` moves only to a
  row inside the table). Loading a file that is not a form is the defect the two last scenarios
  describe (core#1138): a file that is not JSON throws, and a JSON file that is not a form (a saved
  table layout) empties the form. The files are fixtures the project
  already keeps (`fixtures/cars-small.csv`, `fixtures/nx/spgi-legend-position-old.layout`).

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

  @edge-cases @known-failure:39 @realizes:core#1138
  Scenario: A file that is not JSON loaded as a form leaves the form as it was and logs no error
    core#1138: the form decodes the file with no guard, so a CSV file throws "is not valid JSON" and
    the user sees no message. Remove the tag when core#1138 is fixed and such a file is refused.
    When user uploads "fixtures/cars-small.csv" through folder open icon in form viewer
    Then the "fields shown" reading of form viewer should be 11
    And no errors should have been logged

  @edge-cases @known-failure:47 @realizes:core#1138
  Scenario: A JSON file that is not a form loaded as a form leaves the form as it was
    core#1138: the form removes every field before it reads the file, so a JSON file with no form in
    it (here a saved table layout) leaves an empty form with no message. Remove the tag when
    core#1138 is fixed and the fields stay.
    When user uploads "fixtures/nx/spgi-legend-position-old.layout" through folder open icon in form viewer
    Then the "fields shown" reading of form viewer should be 11
