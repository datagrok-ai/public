@viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:viewers.histogram @realizes:viewers.bar-chart @realizes:viewers.pie-chart @realizes:viewers.box-plot @realizes:viewers.pc-plot
Feature: Row Source on a viewer rebound to another table
  The md's second half: each of the seven viewers, added to demog-1000's view, is set to the
  spgi-100 table from its own properties with the md's spgi columns and the
  `${Stereo Category} in ["R_ONE", "S_UNKN"]` filter, and the row sources then have to answer
  spgi-100's filter, selection, current row and hovered group instead of demog-1000's.
  spgi-100: 100 rows, 54 of them R_ONE (36) or S_UNKN (18); row 1 is R_ONE, row 10 is S_ACHIR.
  The md's Filter Panel filter is a TPSA range card on spgi-100 (30 to 80 keeps 59 rows, 32 of the
  54); the selection spans two categories, R_ONE and S_ACHIR (70 rows, 36 of them R_ONE, 22 of
  those inside the TPSA range). The card, the selection and the current row belong to spgi-100 and
  are set on it from demog-1000's view, where the viewer stays on screen (the card goes to the
  filter panel of spgi-100's own view). A one-column viewer repeats its one column pair. With the
  selection cleared, SelectedOrCurrent falls back to the current row — row 1, R_ONE, which the
  filter keeps — and draws nothing once row 10, S_ACHIR, is current. The viewer is a fresh one set
  to spgi-100 through its properties (the md rebinds the viewer that ran the demog half, through
  Context Panel > Data).
  The card stays on to the end: the hovered-group source — a pie chart on Stereo Category (a bar
  chart split by it for the pie chart), bound to spgi-100 the same way — is on Filtered, so its
  R_ONE group is the 22 R_ONE rows inside the TPSA range.
  Not translated: MouseOverRow on the rebound viewer — its hover source is spgi-100's grid, which
  is in a view that is not shown while the viewer is.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens spgi dataset
    And user switches to the "demog-1000" table view

  Scenario Outline: <viewer> rebound to spgi-100 answers that table's rows under every row source
    Given user adds a <viewer> viewer
    When user sets properties of <viewer> viewer:
      | table        | spgi-100                                  |
      | <property>   | <column>                                  |
      | <property 2> | <column 2>                                |
      | rowSource    | Filtered                                  |
      | filter       | ${Stereo Category} in ["R_ONE", "S_UNKN"] |
    Then <viewer> viewer should be bound to table "spgi-100"
    And <viewer> viewer should show 54 rows
    When user adds a range filter on "TPSA" of table "spgi-100" from 30 to 80
    Then 59 rows of table "spgi-100" should pass the filter
    And <viewer> viewer should show 32 rows
    When user sets "rowSource" property of <viewer> viewer to "All"
    Then <viewer> viewer should show 54 rows
    When user sets "rowSource" property of <viewer> viewer to "Selected"
    Then <viewer> viewer should show 0 rows
    When user sets "rowSource" property of <viewer> viewer to "FilteredSelected"
    Then <viewer> viewer should show 0 rows
    When user sets "rowSource" property of <viewer> viewer to "Selected"
    And user selects rows of table "spgi-100" where "Stereo Category" is one of "R_ONE, S_ACHIR"
    Then 70 rows of table "spgi-100" should be selected
    And <viewer> viewer should show 36 rows
    When user sets "rowSource" property of <viewer> viewer to "FilteredSelected"
    Then <viewer> viewer should show 22 rows
    When user sets "rowSource" property of <viewer> viewer to "SelectedOrCurrent"
    Then <viewer> viewer should show 36 rows
    When user clears the row selection of table "spgi-100"
    Then <viewer> viewer should show 1 rows
    When user makes row 10 of table "spgi-100" current
    Then "Stereo Category" of the current row of table "spgi-100" should be "S_ACHIR"
    And <viewer> viewer should show 0 rows
    When user sets "rowSource" property of <viewer> viewer to "CurrentRow"
    Then <viewer> viewer should show 0 rows
    When user makes row 1 of table "spgi-100" current
    Then "Stereo Category" of the current row of table "spgi-100" should be "R_ONE"
    And <viewer> viewer should show 1 rows
    Given user adds a <source> viewer with:
      | table             | spgi-100        |
      | <source property> | Stereo Category |
    When user sets "rowSource" property of <viewer> viewer to "MouseOverGroup"
    And user moves the pointer away from <source> viewer
    Then <viewer> viewer should show 0 rows
    When user hovers over the "<source area> R_ONE" area of <source> viewer
    Then <viewer> viewer should show 22 rows
    When user hovers over the "<source area> S_ACHIR" area of <source> viewer
    Then <viewer> viewer should show 0 rows
    And no errors should have been logged

    Examples:
      | viewer       | property            | column                                   | property 2         | column 2                                 | source    | source property    | source area |
      | scatter plot | xColumnName         | Chemical Space X                         | yColumnName        | Chemical Space Y                         | pie chart | categoryColumnName | slice       |
      | line chart   | xColumnName         | Chemical Space X                         | yColumnNames       | TPSA                                     | pie chart | categoryColumnName | slice       |
      | histogram    | valueColumnName     | TPSA                                     | valueColumnName    | TPSA                                     | pie chart | categoryColumnName | slice       |
      | bar chart    | valueColumnName     | TPSA                                     | splitColumnName    | Stereo Category                          | pie chart | categoryColumnName | slice       |
      | pie chart    | categoryColumnName  | Stereo Category                          | categoryColumnName | Stereo Category                          | bar chart | splitColumnName    | bar         |
      | box plot     | categoryColumnNames | Stereo Category                          | valueColumnName    | TPSA                                     | pie chart | categoryColumnName | slice       |
      | pc plot      | columnNames         | Chemical Space X, Chemical Space Y, TPSA | columnNames        | Chemical Space X, Chemical Space Y, TPSA | pie chart | categoryColumnName | slice       |
