@journey @realizes:chem.cp.sar-matrix
Feature: SAR Matrix Summary on series cut from the molecules
  The SAR Matrix demo project: 134 series cut from 10,000 molecules with hERG pIC50, in two fold levels,
  restored as saved rather than rebuilt. Every scenario starts by putting its own tab on screen.
  The grid icon on a key series opens the series on the SAR Matrix tab with its best measured compound
  current and in the Context Panel; on an analog it opens the analog's series with its prediction in
  the Context Panel. A level chip reads the whole Summary at that level and All reads it over every
  series again. The cart and Add all collect analogs into the Make list without doubles; the Make list
  removes a picked compound, hands a copy of itself to the workspace and empties on Clear.
  On the SAR Matrix tab a measured cell makes its compound current and opens its SAR analysis, a
  predicted one opens the prediction; Ctrl-click and Shift-click add the compounds of measured cells to
  the table's selection; the matrix's cart, the Add to Make list button of a cell's SAR analysis and the
  cell menu each add a predicted analog once; Caption prints
  the substituents' molecular weight; Rank by reorders the series list without a rebuild and the series
  filter counts what it keeps. On the SAR Transfer tab a filter that hides the transfer on screen moves the
  pane to the one it keeps, and one that hides every transfer says so.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the context panel is open
    And user opens the "SarMatrixDemo" project
    And SAR Matrix Viewer viewer should have built its matrices

  Scenario: A key series opens its best measured compound on the SAR Matrix tab
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    Then the "restored" reading of SAR Matrix Viewer viewer should be "true"
    And the "matrices" reading of SAR Matrix Viewer viewer should be 134
    When user clicks on open icon of first key series row
    Then the "tab" reading of SAR Matrix Viewer viewer should be "SAR Matrix"
    And "CMPD_CHEMBLID" of the current row should be "CHEMBL3422978"
    And context panel should contain text "SAR analysis"
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    Then the "series" reading of SAR Matrix Viewer viewer should be the text of name of first key series row
    And no errors should have been logged

  Scenario: An analog opens its prediction on the SAR Matrix tab
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    And user clicks on open icon of first analog worth making row
    Then the "tab" reading of SAR Matrix Viewer viewer should be "SAR Matrix"
    And context panel should contain text "Virtual analog"
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    Then name of first analog worth making row should contain the "series" reading of SAR Matrix Viewer viewer
    And no errors should have been logged

  Scenario: A level chip reads the Summary at one level
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    Then All level chip should be selected
    When user remembers the "key series" reading of SAR Matrix Viewer viewer
    And user clicks on L2 level chip
    Then L2 level chip should be selected
    And All level chip should not be selected
    And the "key series" reading of SAR Matrix Viewer viewer should not be as remembered
    When user clicks on All level chip
    Then All level chip should be selected
    And the "key series" reading of SAR Matrix Viewer viewer should be as remembered
    And no errors should have been logged

  Scenario: The Make list collects without doubles, removes, hands a copy to the workspace and clears
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    And user clicks on cart of first analog worth making row
    Then an info balloon containing "Added 1 to the Make list (1 total)" should have been shown
    When user clicks on cart of first analog worth making row
    Then an info balloon containing "Added 0 to the Make list, 1 already there (1 total)" should have been shown
    And Add all to Make list button should have text "Add all 14 to Make list"
    When user clicks on Add all to Make list button
    Then an info balloon containing "Added 13 to the Make list, 1 already there (14 total)" should have been shown
    When user clicks on Make list tab
    Then Make list title should have text "14 compounds"
    When user clicks on the "cell 1 of Status" area of Make list grid viewer
    And user clicks on Remove button
    Then Make list title should have text "13 compounds"
    When user clicks on Add to workspace button
    Then the "SAR make list" table view should open with 13 rows
    When user switches to the "SAR Matrix" view
    And user clicks on Make list tab
    And user clicks on Clear button
    Then SAR Matrix note should contain text "Nothing collected yet"
    And the "make list" reading of SAR Matrix Viewer viewer should be 0
    And no errors should have been logged

  Scenario: A measured cell opens its compound and a predicted cell its prediction
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    And user clicks on first series card
    Then the "series" reading of SAR Matrix Viewer viewer should be the text of name of first series card
    When user clicks on the "cell 3 of c0" area of grid viewer
    Then "CMPD_CHEMBLID" of the current row should be "CHEMBL3422763"
    And context panel should contain text "SAR analysis"
    And context panel should not contain text "Virtual analog"
    When user clicks on the "cell 1 of c0" area of grid viewer
    Then context panel should contain text "Virtual analog"
    And no errors should have been logged

  Scenario: Ctrl-click and Shift-click add measured compounds to the selection
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    And user clicks on first series card
    And user clears the row selection
    And user clicks on the "cell 3 of c0" area of grid viewer holding Control
    Then 1 row should be selected
    When user clicks on the "cell 4 of c0" area of grid viewer holding Control
    Then 2 rows should be selected
    When user clicks on the "cell 5 of c0" area of grid viewer holding Shift
    Then 3 rows should be selected
    When user clicks on the "cell 4 of c0" area of grid viewer holding Control
    Then 2 rows should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: The matrix's cart, a cell's SAR analysis and the cell menu each add a predicted analog once
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    And user clicks on first series card
    And user clicks on the "cell 1 of c0" area of grid viewer
    And user clicks on matrix cart
    Then an info balloon containing "Virtual analog added to the Make list (1 total)" should have been shown
    And the "make list" reading of SAR Matrix Viewer viewer should be 1
    When user clicks on the "cell 6 of c0" area of grid viewer
    And user clicks on Add to Make list button in context panel
    Then an info balloon containing "Virtual analog added to the Make list (2 total)" should have been shown
    And the "make list" reading of SAR Matrix Viewer viewer should be 2
    When user picks "SAR Matrix > Add this compound to make list" from the context menu of the "cell 9 of c0" area of grid viewer
    Then the "make list" reading of SAR Matrix Viewer viewer should be 3
    When user picks "SAR Matrix > Add this compound to make list" from the context menu of the "cell 1 of c0" area of grid viewer
    Then an info balloon containing "Virtual analog is already on the Make list (3 total)" should have been shown
    And the "make list" reading of SAR Matrix Viewer viewer should be 3
    When user clicks on Make list tab
    And user clicks on Clear button
    Then the "make list" reading of SAR Matrix Viewer viewer should be 0
    And no errors should have been logged

  Scenario: Caption prints the substituents' molecular weight under the columns
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    And user clicks on first series card
    Then the "column captions" reading of SAR Matrix Viewer viewer should be "μ 7.2, μ 5.6, μ 6.9, μ 7.5, μ 7.3, μ 7.9, μ 8.2, μ 5.0, μ 5.0, μ 5.2, μ 5.1, μ 7.7, μ 7.7, μ 6.1, μ 6.8, μ 5.7, μ 4.9, μ 5.2, μ 4.6, μ 4.6"
    When user selects "Molecular weight" in "Caption" input
    Then the "column captions" reading of SAR Matrix Viewer viewer should be "MW 32, MW 17, MW 16, MW 20, MW 27, MW 36, MW 70, MW 31, MW 32, MW 18, MW 74, MW 46, MW 30, MW 60, MW 80, MW 52, MW 45, MW 46, MW 62, MW 46"
    When user selects "Potency" in "Caption" input
    Then the "column captions" reading of SAR Matrix Viewer viewer should be "μ 7.2, μ 5.6, μ 6.9, μ 7.5, μ 7.3, μ 7.9, μ 8.2, μ 5.0, μ 5.0, μ 5.2, μ 5.1, μ 7.7, μ 7.7, μ 6.1, μ 6.8, μ 5.7, μ 4.9, μ 5.2, μ 4.6, μ 4.6"
    And no errors should have been logged

  Scenario: Rank by reorders the series list without a rebuild, and the series filter counts what it keeps
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    And user remembers the "series list" reading of SAR Matrix Viewer viewer
    And user selects "SAR discontinuity" in "Rank by" input
    Then the "series list" reading of SAR Matrix Viewer viewer should not be as remembered
    And the "restored" reading of SAR Matrix Viewer viewer should be "true"
    When user selects "Potent compounds" in "Rank by" input
    Then the "series list" reading of SAR Matrix Viewer viewer should be as remembered
    When user clicks on series filter icon
    And user hovers over "Series" filter card
    And user clicks on search icon of "Series" filter card
    And user types "Series 69" into the search of the "Series" filter card
    Then series match count should have text "1 of 134 match"
    And the "series list" reading of SAR Matrix Viewer viewer should be "Series 69"
    When user clears the search of the "Series" filter card
    Then the "series list" reading of SAR Matrix Viewer viewer should be as remembered
    And series match count should have text ""
    When user presses Escape
    Then no errors should have been logged

  Scenario: A transfer filter moves the pane to the transfer it keeps, and says when it keeps none
    When user clicks on SAR Transfer tab in SAR Matrix Viewer viewer
    Then transfer title should have text "Series 69 · Core 3 → Series 33 · Core 1"
    When user clicks on transfer filter icon
    And user hovers over "Target" filter card
    And user clicks on search icon of "Target" filter card
    And user types "Series 98" into the search of the "Target" filter card
    Then transfer title should have text "Series 69 · Core 3 → Series 98 · Core 1"
    When user types "Series 5" into the search of the "Target" filter card
    Then SAR Matrix note should have text "No transfer matches the filter — clear it to see the transfers again."
    When user clears the search of the "Target" filter card
    And user presses Escape
    And user selects "Series 33 · Core 1 · R1" in "Transfer to" input
    Then transfer title should have text "Series 69 · Core 3 → Series 33 · Core 1"
    And no errors should have been logged
