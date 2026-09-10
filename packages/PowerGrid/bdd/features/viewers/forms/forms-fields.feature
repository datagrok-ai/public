@journey @viewers @realizes:viewers.forms
Feature: Forms viewer field set and number format
  The life of the field set: the Fields property picks the columns and their order, the header's
  remove icon drops one, a column renamed to a "~" name leaves the set, an empty set draws nothing
  and says nothing, a column removed from the table takes its field with it, and a named number
  format reaches the float fields and leaves the integer and string ones alone. One journey on
  demog-1000 with Show Mouse Over Row off.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer with:
      | Show Mouse Over Row | false |
    Then forms viewer should be visible
    And the "fields shown" reading of forms viewer should be 11

  Scenario: The fields are the picked columns, in the picked order
    When user sets "Fields" property of forms viewer to "RACE, AGE, SEX"
    Then the "fields" reading of forms viewer should be "RACE, AGE, SEX"
    And the "fields shown" reading of forms viewer should be 3
    And "Fields" property of forms viewer should be "RACE, AGE, SEX"
    And forms viewer should have a "label RACE" area
    And forms viewer should have a "label AGE" area
    And no errors should have been logged

  Scenario: The header's remove icon drops a field
    When user clicks on the "remove AGE" area of forms viewer
    Then the "fields" reading of forms viewer should be "RACE, SEX"
    And the "fields shown" reading of forms viewer should be 2
    And forms viewer should not have a "label AGE" area
    And no errors should have been logged

  Scenario: A column renamed to a "~" name leaves the field set
    When user renames "RACE" column to "~RACE"
    Then the "fields" reading of forms viewer should be "SEX"
    And the "fields shown" reading of forms viewer should be 1
    And forms viewer should not have a "label RACE" area
    When user renames "~RACE" column to "RACE"
    Then the table should have a column "RACE"
    And no errors should have been logged

  Scenario: A named number format reaches the float fields only
    When user sets "Fields" property of forms viewer to "HEIGHT, AGE, SEX"
    And user makes row 1 current
    Then "Number Format" property of forms viewer should be "Same as grid"
    And the "AGE of current card" reading of forms viewer should be "26"
    And the "SEX of current card" reading of forms viewer should be "F"
    When user sets "Number Format" property of forms viewer to "3 significant digits"
    Then the "HEIGHT of current card" reading of forms viewer should be "175"
    And the "AGE of current card" reading of forms viewer should be "26"
    And the "SEX of current card" reading of forms viewer should be "F"
    When user sets "Number Format" property of forms viewer to "3 digits after comma"
    Then the "HEIGHT of current card" reading of forms viewer should be "174.705"
    When user sets "Number Format" property of forms viewer to "Same as grid"
    Then "Number Format" property of forms viewer should be "Same as grid"
    And no errors should have been logged

  Scenario: An empty field set draws no card and says nothing
    When user sets "Fields" property of forms viewer to ""
    Then the "fields shown" reading of forms viewer should be 0
    And the "cards" reading of forms viewer should be 0
    And forms viewer should not have a "label HEIGHT" area
    And no error or warning balloon should have been shown
    When user sets "Fields" property of forms viewer to "USUBJID, AGE, SEX"
    Then the "fields shown" reading of forms viewer should be 3
    And no errors should have been logged

  Scenario: A column removed from the table takes its field with it
    When user removes "AGE" column
    Then the table should not have a column "AGE"
    And the "fields" reading of forms viewer should be "USUBJID, SEX"
    And forms viewer should not have a "label AGE" area
    And no error or warning balloon should have been shown
    And no errors should have been logged
