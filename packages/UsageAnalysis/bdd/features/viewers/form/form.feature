@journey @viewers @realizes:viewers.form
Feature: Form viewer
  One record as a form: every column of the current row as a labelled field, the toolbar that walks
  the rows and selects them, the track-row modes, editing a field, the field set and its
  persistence, and the colour coding the grid gives a field. One journey on demog-1000; every
  scenario puts back what it changed. The viewer draws no canvas, so its claims are the readings
  and the hit areas it reports, never pixels. A default form orders its fields by a relevance
  score (a molecule first, then anything with a semantic type, constant columns last) and
  re-sorts a set it is given the same way, so the claims name fields and count them, never
  their order. A column removed from the table stays on a designed form with an empty field.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a form viewer
    Then form viewer should be visible
    And the "fields shown" reading of form viewer should be 11
    And form viewer should have a "field USUBJID" area
    And form viewer should have a "field SEVERITY" area

  Scenario: The fields show the current row
    Then the "row" reading of form viewer should be 1
    And the "USUBJID" reading of form viewer should be "X0273T21000300003"
    And the "AGE" reading of form viewer should be "26"
    And the "SEX" reading of form viewer should be "F"
    And the "RACE" reading of form viewer should be "Caucasian"
    And the "DIS_POP" reading of form viewer should be "Indigestion"
    And form viewer should have a "field AGE" area
    And form viewer should have a "label AGE" area
    And no errors should have been logged

  Scenario: A column's colour coding reaches its field
    Then the "background of AGE" reading of form viewer should be ""
    When user colors "AGE" column linearly from "#FF0000" to "#0000FF"
    Then the "background of AGE" reading of form viewer should differ from before
    When user removes the coloring of "AGE" column
    Then the "background of AGE" reading of form viewer should be ""
    And no errors should have been logged

  Scenario: The arrows walk the rows
    When user clicks on the "next row" area of form viewer
    Then row 2 should be current
    And the "USUBJID" reading of form viewer should be "X0273T21000300005"
    And the "AGE" reading of form viewer should be "30"
    When user clicks on the "previous row" area of form viewer
    Then row 1 should be current
    And the "USUBJID" reading of form viewer should be "X0273T21000300003"
    And no errors should have been logged

  Scenario: The form follows the grid's current row
    When user makes row 43 current
    Then the "row" reading of form viewer should be 43
    And the "USUBJID" reading of form viewer should be "X0273T21002400007"
    And the "AGE" reading of form viewer should be "54"
    When user makes row 1 current
    Then the "row" reading of form viewer should be 1

  Scenario: The row selector toggles the row's selection
    Given user clears the row selection
    When user clicks on the "row selector" area of form viewer
    Then only rows where "USUBJID" is "X0273T21000300003" should be selected
    And the "row selected" reading of form viewer should be "true"
    And "square" icon in form viewer should be selected
    When user clicks on the "row selector" area of form viewer
    Then no rows should be selected
    And the "row selected" reading of form viewer should be "false"
    And "square" icon in form viewer should not be selected

  Scenario: Track Row: the form follows the mouse-over row, or nothing
    When user picks "Track Row > Mouse Over" from the context menu of form viewer
    Then "Sync Mode" property of form viewer should be "Mouse Over"
    When user hovers over the "cell 2 of AGE" area of grid
    Then the "USUBJID" reading of form viewer should be "X0273T21000300005"
    When user moves the pointer away from grid
    And user picks "Track Row > None" from the context menu of form viewer
    Then "Sync Mode" property of form viewer should be "None"
    When user makes row 5 current
    Then the "USUBJID" reading of form viewer should be "X0273T21000300005"
    When user picks "Track Row > Current" from the context menu of form viewer
    Then "Sync Mode" property of form viewer should be "Current"
    And the "USUBJID" reading of form viewer should be "X0273T21000500006"
    When user makes row 1 current

  Scenario: The keyboard walks and selects rows
    Given user clears the row selection
    When user focuses on form viewer
    And user presses ArrowRight in form viewer
    Then row 2 should be current
    When user presses ArrowDown in form viewer
    Then row 3 should be current
    When user presses ArrowLeft in form viewer
    Then row 2 should be current
    When user presses ArrowUp in form viewer
    Then row 1 should be current
    When user presses Space in form viewer
    Then only rows where "USUBJID" is "X0273T21000300003" should be selected
    When user clears the row selection

  Scenario: A field writes its cell only while the form is editable
    Then the "editable" reading of form viewer should be "false"
    When user clicks on the "edit" area of form viewer
    Then the "editable" reading of form viewer should be "true"
    When user enters "31" into the "field AGE" area of form viewer
    Then the value of "AGE" column in row 1 should be "31"
    And the "AGE" reading of form viewer should be "31"
    When user enters "26" into the "field AGE" area of form viewer
    And user clicks on the "edit" area of form viewer
    Then the "editable" reading of form viewer should be "false"
    And the value of "AGE" column in row 1 should be "26"
    When user enters "99" into the "field AGE" area of form viewer
    Then the value of "AGE" column in row 1 should be "26"

  Scenario: The toolbar shows what the properties allow
    When user sets "Show Next Row Arrow" property of form viewer to "false"
    Then form viewer should not have a "next row" area
    And form viewer should have a "previous row" area
    When user sets "Show Row Selector" property of form viewer to "false"
    Then form viewer should not have a "row selector" area
    When user sets properties of form viewer:
      | Show Next Row Arrow | true |
      | Show Row Selector   | true |
    Then form viewer should have a "next row" area
    And form viewer should have a "row selector" area
    When user sets "Show Navigation" property of form viewer to "false"
    Then form viewer should not have a "next row" area
    And form viewer should not have a "select columns" area
    When user sets "Show Navigation" property of form viewer to "true"
    Then form viewer should have a "next row" area

  Scenario: The field set follows the columns it is given
    When user sets "columnNames" property of form viewer to "AGE, SEX"
    Then the "fields shown" reading of form viewer should be 2
    And form viewer should have a "field SEX" area
    And form viewer should have a "field AGE" area
    And form viewer should not have a "field RACE" area
    And the "AGE" reading of form viewer should be "26"

  Scenario: The field set survives a layout round-trip
    Given user sets "columnNames" property of form viewer to "AGE, HEIGHT, WEIGHT"
    And the "fields shown" reading of form viewer should be 3
    When user saves the layout of the current table view
    And user clicks on close icon of form viewer
    Then form viewer should be absent
    When user loads the saved layout
    Then form viewer should be visible
    And the "fields shown" reading of form viewer should be 3
    And form viewer should have a "field HEIGHT" area
    And form viewer should have a "field WEIGHT" area
    And the "AGE" reading of form viewer should be "26"
    And no errors should have been logged

  Scenario: A removed column empties its field but keeps it on the form
    When user removes "WEIGHT" column
    Then the table should not have a column "WEIGHT"
    And the "fields shown" reading of form viewer should be 3
    And the "WEIGHT" reading of form viewer should be ""
    And the "AGE" reading of form viewer should be "26"
    And no errors should have been logged
    When user adds a calculated column "WEIGHT" with formula "0"
    Then the table should have a column "WEIGHT"

  Scenario: The arrows walk the filtered rows
    Given user sets "columnNames" property of form viewer to "USUBJID, AGE, SEX"
    When user filters rows where "SEX" is "M"
    And user makes row 4 current
    Then the "USUBJID" reading of form viewer should be "X0273T21000400002"
    When user clicks on the "next row" area of form viewer
    Then the "SEX" reading of form viewer should be "M"
    And the "USUBJID" reading of form viewer should be "X0273T21000500008"
    When user resets the filter
    And user makes row 1 current

  Scenario: The viewer closes from its title bar
    When user clicks on close icon of form viewer
    Then form viewer should be absent
    And no errors should have been logged
