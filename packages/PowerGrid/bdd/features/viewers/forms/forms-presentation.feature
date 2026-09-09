@journey @viewers @realizes:viewers.forms
Feature: Forms viewer presentation
  What a field looks like: a text column is an input, a column's colour coding paints the field's
  background on every card that shows the column, Color Code off takes the paint away and back on
  restores exactly the colour the scheme gives, and removing the coloring clears it. One journey on
  demog-1000 with Show Mouse Over Row off. The renderer-size ladder and the molecule and curve
  fields need Chem and Curves on the stand and are not translated yet.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer with:
      | Show Mouse Over Row | false |
    And user makes row 1 current
    Then the "USUBJID of current card" reading of forms viewer should be "X0273T21000300003"

  Scenario: A text column is an input field
    Then the "field kind of USUBJID" reading of forms viewer should be "input"
    And the "field kind of AGE" reading of forms viewer should be "input"
    And no errors should have been logged

  Scenario: A column's colour coding paints its field and Color Code gates it
    When user colors "AGE" column linearly from "#FF0000" to "#0000FF"
    Then the "background of AGE of current card" reading of forms viewer should differ from before
    When user remembers the "background of AGE of current card" reading of forms viewer
    And user sets "Color Code" property of forms viewer to "false"
    Then the "background of AGE of current card" reading of forms viewer should differ from before
    When user sets "Color Code" property of forms viewer to "true"
    Then the "background of AGE of current card" reading of forms viewer should be as remembered
    When user removes the coloring of "AGE" column
    Then the "background of AGE of current card" reading of forms viewer should differ from before
    And no errors should have been logged

  Scenario: The colouring reaches the selected rows' cards too
    When user selects rows where "USUBJID" is one of "X0273T21000400001, X0273T21001500015"
    Then the "cards" reading of forms viewer should be 3
    When user colors "AGE" column linearly from "#FF0000" to "#0000FF"
    Then the "background of AGE of card 2" reading of forms viewer should differ from before
    And the "background of AGE of card 3" reading of forms viewer should differ from before
    When user removes the coloring of "AGE" column
    And user clears the row selection
    Then the "cards" reading of forms viewer should be 1
    And no errors should have been logged
