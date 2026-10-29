@journey @viewers @realizes:viewers.form
Feature: Form viewer
  One record as a form: every column of the current row as a labelled field, the toolbar that walks
  the rows and selects them, the track-row modes, editing a field, the field set and its
  persistence (a layout, and a form saved to a file and loaded back from it through the toolbar's
  own icons), and the colour coding the grid gives a field. One journey on demog-1000; every
  scenario puts back what it changed. The viewer draws no canvas, so its claims are the readings
  and the hit areas it reports, never pixels. A default form orders its fields by a relevance
  score (a molecule first, then anything with a semantic type, constant columns last) and
  re-sorts a set it is given the same way, so the claims name fields and count them, never
  their order.
  The md's manual design mode is a scenario of its own: with the toolbar's design toggle on, a field
  dragged down stays below the fields it passed once the toggle is off, and a form saved to a file
  after that brings the field back to its new place when the file is loaded over a rebuilt form.
  From the md's manual sections: Column selector — the toolbar's column list dialog builds the form
  from the columns checked in it, and opened again it starts from the fields shown; Column changes
  reaction — a column renamed in the grid's Column Properties dialog renames its field and the
  label beside it (step 3, a deleted column, is claimed in `form-edge-cases`);
  Table switching — the Table property moves a default form to spgi-100's columns and back, and a
  form built from columns chosen in that dialog shows spgi-100's columns too.

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

  Scenario: Bound to another table the form shows that table's columns
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "table" property of form viewer to "spgi-100"
    Then form viewer should be bound to table "spgi-100"
    And the "fields shown" reading of form viewer should be 13
    And form viewer should have a "field TPSA" area
    And form viewer should have a "field Series" area
    And form viewer should not have a "field AGE" area
    When user sets "table" property of form viewer to "demog-1000"
    Then form viewer should be bound to table "demog-1000"
    And the "fields shown" reading of form viewer should be 11
    And form viewer should have a "field USUBJID" area
    And form viewer should not have a "field TPSA" area
    And the "AGE" reading of form viewer should be "26"
    And no errors should have been logged

  Scenario: A form built from chosen columns shows the new table's columns once Table is switched
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user clicks on the "select columns" area of form viewer
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    And user toggles the "AGE" column in the column list of "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "fields shown" reading of form viewer should be 2
    When user sets "table" property of form viewer to "spgi-100"
    Then form viewer should be bound to table "spgi-100"
    And the "fields shown" reading of form viewer should be 13
    And form viewer should have a "field TPSA" area
    And form viewer should not have a "field AGE" area
    When user sets "table" property of form viewer to "demog-1000"
    Then form viewer should be bound to table "demog-1000"
    And the "fields shown" reading of form viewer should be 11
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

  Scenario: The toolbar's column list builds the form from the columns checked in it
    When user clicks on the "select columns" area of form viewer
    Then "Select columns..." dialog should be visible
    And "11 checked" text in "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user toggles the "AGE" column in the column list of "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    And the "fields shown" reading of form viewer should be 2
    And form viewer should have a "field AGE" area
    And form viewer should have a "field SEX" area
    And form viewer should not have a "field RACE" area
    And the "AGE" reading of form viewer should be "26"
    When user clicks on the "select columns" area of form viewer
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user toggles the "RACE" column in the column list of "Select columns..." dialog
    Then "3 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    And the "fields shown" reading of form viewer should be 3
    And form viewer should have a "field RACE" area
    And form viewer should have a "field SEX" area
    And the "RACE" reading of form viewer should be "Caucasian"
    And no errors should have been logged

  Scenario: A column renamed in the grid renames its field and the field's label
    Then form viewer should have a "label AGE" area
    When user picks "Column Properties..." from the context menu of the "header AGE" area of grid
    Then "AGE" dialog should be visible
    When user types "AGE_NEW" into "New name:" input in "AGE" dialog
    And user clicks on OK button in "AGE" dialog
    Then the table should have a column "AGE_NEW"
    And form viewer should have a "field AGE_NEW" area
    And form viewer should have a "label AGE_NEW" area
    And form viewer should not have a "label AGE" area
    And form viewer should not have a "field AGE" area
    And the "AGE_NEW" reading of form viewer should be "26"
    When user renames "AGE_NEW" column to "AGE"
    Then form viewer should have a "field AGE" area
    And form viewer should have a "label AGE" area
    And no errors should have been logged

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

  Scenario: A form saved to a file comes back from the file
    Given user sets "columnNames" property of form viewer to "AGE, SEX, RACE"
    Then the "fields shown" reading of form viewer should be 3
    When user downloads a file through arrow to bottom icon in form viewer
    Then the downloaded file should contain "RACE"
    And the downloaded file should not contain "HEIGHT"
    When user sets "columnNames" property of form viewer to "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    Then the "fields shown" reading of form viewer should be 11
    When user uploads the downloaded file through folder open icon in form viewer
    Then the "fields shown" reading of form viewer should be 3
    And form viewer should have a "field RACE" area
    And form viewer should have a "field SEX" area
    And form viewer should not have a "field HEIGHT" area
    And the "AGE" reading of form viewer should be "26"
    When user sets "columnNames" property of form viewer to "AGE, HEIGHT, WEIGHT"
    Then the "fields shown" reading of form viewer should be 3
    And no errors should have been logged

  Scenario: A field dragged in design mode keeps its place, and the saved file brings the place back
    Given user sets "columnNames" property of form viewer to "AGE, SEX, RACE"
    Then the "fields shown" reading of form viewer should be 3
    And the "design mode" reading of form viewer should be "false"
    And the "field RACE" area of form viewer should hang below the "field AGE" area
    When user clicks on the "design" area of form viewer
    Then the "design mode" reading of form viewer should be "true"
    When user drags the "field AGE" area of form viewer by 160 pixels to the down
    Then the "field AGE" area of form viewer should hang below the "field RACE" area
    And the "field AGE" area of form viewer should hang below the "field SEX" area
    When user clicks on the "design" area of form viewer
    Then the "design mode" reading of form viewer should be "false"
    And the "field AGE" area of form viewer should hang below the "field RACE" area
    When user downloads a file through arrow to bottom icon in form viewer
    And user sets "columnNames" property of form viewer to "AGE, SEX, RACE, HEIGHT"
    Then the "fields shown" reading of form viewer should be 4
    And the "field RACE" area of form viewer should hang below the "field AGE" area
    When user uploads the downloaded file through folder open icon in form viewer
    Then the "fields shown" reading of form viewer should be 3
    And the "field AGE" area of form viewer should hang below the "field RACE" area
    When user sets "columnNames" property of form viewer to "AGE, HEIGHT, WEIGHT"
    Then the "fields shown" reading of form viewer should be 3
    And no errors should have been logged

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

