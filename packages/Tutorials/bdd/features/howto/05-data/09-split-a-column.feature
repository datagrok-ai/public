@guide @tutorial
Feature: How-to · Data
  For the maintainer: splitting a text column two ways, in one clip (it took in the former clip "Split a column by a
  pattern"). demog's DEMOG column packs three values into one text ("53 C F": age, the first letter of race, sex).
  With a DEMOG cell current, Data > Split Column opens on that column; a space as the separator shows col1, col2 and
  col3 in the preview, and OK adds them to the table (col1 stays text: one row of demog reads ". C M", its age
  missing as "."). Then the titanic tab: its name column reads "Surname, Title. Given names" ("Allen, Miss.
  Elisabeth Walton"), where no one separator does, and Data > Split Column by RegExp with the pattern
  (.+), (\w+)\. (.+) adds one column per group in parentheses: r the surname, r (2) the title, r (3) the given
  names, previewed in a window beside the dialog while the pattern is typed. Both tables are opened off camera,
  titanic first (its grid shows five columns, set off camera so the new ones land in view), so demog is the view
  the clip starts on and titanic's tab is beside it; nothing is made on the server.

  Scenario: Split a column by a separator or a pattern
    Split a text column into several at a separator, or pull its parts out with a pattern.

    Given user is logged in
    And simple mode is off
    And user opens titanic dataset
    And the grid shows only the columns "pclass, survived, name, sex, age"
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "DEMOG holds three values in one text: age, race, sex"
    # caption: Click the column to split
    And user clicks on the "cell 1 of DEMOG" area of grid
    # caption: Split Column: one text column into several
    And user picks "Data > Split Column..." from the top menu
    Then Split Column dialog should be visible
    And Column input in Split Column dialog should contain text "DEMOG"
    # caption: Split by: here a space
    When user enters " " into "Split by" input in Split Column dialog
    And user clicks on OK button in Split Column dialog
    Then Split Column dialog should be hidden
    And the table should have a column "col3"
    And the value of "col1" column in row 1 should be "53"
    And the value of "col2" column in row 1 should be "C"
    And the value of "col3" column in row 1 should be "F"
    When the video says "Three new columns: age, race and sex" at the bottom
    And user hovers over the "header col1" area of grid
    When user clicks on the tab of the "titanic" view
    # caption: In titanic, each name holds a surname, a title and given names
    And user hovers over the "cell 1 of name" area of grid
    # caption: No one separator splits these: split by a pattern
    And user picks "Data > Split Column by RegExp..." from the top menu
    Then Split Column by RegExp dialog should be visible
    When user selects "name" in "Column" input in Split Column by RegExp dialog
    # caption: Each (group) in the pattern becomes a column
    And user enters "(.+), (\w+)\. (.+)" into "RegExp" input in Split Column by RegExp dialog
    # caption: The preview shows the parts it found
    And user hovers over "Preview result columns" dialog
    And user clicks on OK button in Split Column by RegExp dialog
    Then Split Column by RegExp dialog should be hidden
    # the code asks for r1, r2, r3 (prefix "r", startWith 1), but ColumnList.getUnusedName ignores startWith
    And the table should have a column "r (3)"
    And the value of "r" column in row 1 should be "Allen"
    And the value of "r (2)" column in row 1 should be "Miss"
    And the value of "r (3)" column in row 1 should be "Elisabeth Walton"
    When the video says "Surname, title and given names, each in a column" at the bottom
    And user hovers over the "cell 1 of r (2)" area of grid
    Then no errors should have been logged
