@journey @viewers @realizes:viewers.forms
Feature: Forms viewer field lifecycle and number format
  The field set is a list of column names the viewer keeps, and this is what happens to it: it is
  drawn in the order it was given, a header ✕ takes one out, a column that leaves the table takes
  its field with it, a rename carries the field to the new name, and a rename to a `~` name — the
  platform's mark for a service column — drops it. Then the number format, which is the viewer's
  own and not the grid's.
  Two readings carry all of it. `fields` is the CONFIGURED set; `header labels` is what the header
  actually drew, and they differ while the viewer has not yet pruned a field whose column is gone —
  which is what makes "a dropped or `~`-renamed column prunes the field" an honest claim rather
  than a repaint.
  The viewer has NO error state. An empty field set draws zero labels and zero fields and shows
  nothing at all: no balloon, no banner, no message. The 20-column cap is silent the same way. So
  the empty-field-set scenario asserts `fields shown` is 0 and `header labels` is empty, and that
  the viewer reports no error — never that it says something.
  The spec this replaces removed a field through a 30-line retry loop around the header ✕ that
  threw when the icon had a zero box, and drove the grid's header context menu through two more
  retry loops to drop and rename a column. The ✕ is the `remove <COL>` hit area, and the column
  steps are the platform's.
  COMPUTED_H is `${HEIGHT}` with no format of its own, so under "Same as grid" the card shows what
  the grid shows (174.71) and under "3 digits after comma" it shows the viewer's own (174.705) —
  the underlying value never moves.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calculated column "COMPUTED_H" with formula "${HEIGHT}"
    And user adds a forms viewer
    Then forms viewer should be visible
    And the table should have 12 columns
    And the "fields shown" reading of forms viewer should be 12
    And the "fields" and "header labels" readings of forms viewer should be the same

  Scenario: The fields are drawn in the order they were given, not in table order
    When user sets "fieldsColumnNames" property of forms viewer to "RACE, AGE, SEX"
    Then the "fields" reading of forms viewer should be "RACE, AGE, SEX"
    And the "header labels" reading of forms viewer should be "RACE, AGE, SEX"
    And the "fields shown" reading of forms viewer should be 3
    And forms viewer should have a "field RACE of card 1" area
    And forms viewer should not have a "field WEIGHT of card 1" area
    And the "RACE of card 1" reading of forms viewer should be "Caucasian"
    When user sets "fieldsColumnNames" property of forms viewer to "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, COMPUTED_H"
    Then the "fields shown" reading of forms viewer should be 12
    And no errors should have been logged

  Scenario: The header cross takes a field out and leaves the order of the rest
    When user sets "fieldsColumnNames" property of forms viewer to "RACE, AGE, SEX"
    Then forms viewer should have a "remove AGE" area
    When user clicks on the "remove AGE" area of forms viewer
    Then the "fields" reading of forms viewer should be "RACE, SEX"
    And the "header labels" reading of forms viewer should be "RACE, SEX"
    And the "fields shown" reading of forms viewer should be 2
    And forms viewer should not have a "field AGE of card 1" area
    And forms viewer should not have a "label AGE" area
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Number Format is the viewer's own and leaves the value alone
    When user sets "fieldsColumnNames" property of forms viewer to "COMPUTED_H, AGE, SEX"
    Then "numberFormat" property of forms viewer should be "Same as grid"
    And the "COMPUTED_H of card 1" reading of forms viewer should be "174.71"
    And the "COMPUTED_H" cell of row 1 should be displayed as "174.71"
    And the "AGE of card 1" reading of forms viewer should be "26"
    And the "SEX of card 1" reading of forms viewer should be "F"
    When user sets "numberFormat" property of forms viewer to "3 digits after comma"
    Then the "COMPUTED_H of card 1" reading of forms viewer should be "174.705"
    And the "COMPUTED_H" cell of row 1 should be displayed as "174.71"
    And the "AGE of card 1" reading of forms viewer should be "26"
    And the "SEX of card 1" reading of forms viewer should be "F"
    And the value of "COMPUTED_H" column in row 1 should be "174.7050018310547"
    When user sets "numberFormat" property of forms viewer to "2 digits after comma"
    Then the "COMPUTED_H of card 1" reading of forms viewer should be "174.71"
    When user sets "numberFormat" property of forms viewer to "Same as grid"
    Then the "COMPUTED_H of card 1" reading of forms viewer should be "174.71"
    And no errors should have been logged

  Scenario: An empty field set draws nothing and says nothing
    When user sets "fieldsColumnNames" property of forms viewer to ""
    Then the "fields shown" reading of forms viewer should be 0
    And the "fields" reading of forms viewer should be ""
    And the "header labels" reading of forms viewer should be ""
    And forms viewer should not have a "label AGE" area
    And forms viewer should report no error
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user sets "fieldsColumnNames" property of forms viewer to "RACE, AGE, SEX"
    Then the "fields shown" reading of forms viewer should be 3

  Scenario: A renamed column carries its field to the new name
    Then the "fields" reading of forms viewer should be "RACE, AGE, SEX"
    When user renames "SEX" column to "GENDER"
    Then the table should have a column "GENDER"
    And the "fields" reading of forms viewer should be "RACE, AGE, GENDER"
    And the "header labels" reading of forms viewer should be "RACE, AGE, GENDER"
    And forms viewer should have a "field GENDER of card 1" area
    And forms viewer should not have a "field SEX of card 1" area
    And the "GENDER of card 1" reading of forms viewer should be "F"
    And no error or warning balloon should have been shown
    When user renames "GENDER" column to "SEX"
    Then the "fields" reading of forms viewer should be "RACE, AGE, SEX"
    And no errors should have been logged

  Scenario: A rename to a service name drops the field
    Then the "fields" reading of forms viewer should be "RACE, AGE, SEX"
    When user renames "SEX" column to "~SERVICE"
    Then the "fields" reading of forms viewer should be "RACE, AGE"
    And the "header labels" reading of forms viewer should be "RACE, AGE"
    And the "header labels" reading of forms viewer should not contain "~"
    And the "fields shown" reading of forms viewer should be 2
    And forms viewer should not have a "field ~SERVICE of card 1" area
    And no error or warning balloon should have been shown
    When user renames "~SERVICE" column to "SEX"
    Then the table should have a column "SEX"
    And no errors should have been logged

  Scenario: A column that leaves the table takes its field with it
    When user sets "fieldsColumnNames" property of forms viewer to "RACE, AGE, SEX"
    Then the "fields" reading of forms viewer should be "RACE, AGE, SEX"
    When user removes "RACE" column
    Then the table should not have a column "RACE"
    And the "fields" reading of forms viewer should be "AGE, SEX"
    And the "header labels" reading of forms viewer should be "AGE, SEX"
    And the "fields shown" reading of forms viewer should be 2
    And forms viewer should not have a "field RACE of card 1" area
    And the "AGE of card 1" reading of forms viewer should be "26"
    And forms viewer should report no error
    And no error or warning balloon should have been shown
    And no errors should have been logged
