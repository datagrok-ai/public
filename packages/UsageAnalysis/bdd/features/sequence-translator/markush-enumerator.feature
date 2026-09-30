@realizes:chem.menu.transform.markush-enumeration
Feature: Markush Enumerator: a core from a cell, R-groups from a table, Cartesian and Zip
  The Markush Enumerator (SequenceTranslator) enumerates molecules from R-labelled cores and
  R-group lists. The core C1C2CN([*:2])CC2CN1[*:1] (row 2 of chem_enum_cores) has two positions,
  R1 and R2; chem_enum_rgroups holds 4 distinct R1 and 4 distinct R2 substituents, so Cartesian
  gives 4 × 4 = 16 molecules and Zip gives 4. The result table takes the name typed into Table name
  (GROK-20223); the Markush Enumerator app keeps its content when the user leaves its view and
  comes back (GROK-20372). Translated from the TestTrack case SequenceTranslator/markush-enumerator.

  The dialog's core and R-group cards have no names or readings: the R-group counts are read from
  the cards' captions ("r group 4 · R1" is there, "r group 5 · R1" is not), and "exactly one core
  card" is kept without (see the request document). Which core the dialog took is proven by the
  result's Core column (every row holds row 2's core; rows 3 and 4 have R1 and R2 too); the top-menu
  scenario ends on CANCEL, so there "that core" is shown only by its R-numbers (R2 present, R3
  absent — true of rows 2, 3 and 4 alike). The two result tables live in the workspace only and go
  with the views at the end of each scenario.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator---tests tree node inside browse tree is expanded
    When user double-clicks Files---App-Data---SequenceTranslator---tests---chem_enum_rgroups.csv tree node inside browse tree
    Then the "chem_enum_rgroups" view should be current
    When user clicks on browse tab
    And user double-clicks Files---App-Data---SequenceTranslator---tests---chem_enum_cores.csv tree node inside browse tree
    Then the "chem_enum_cores" view should be current
    And "Core" column should have semantic type "Molecule"
    And the value of "Core" column in row 2 should be "C1C2CN([*:2])CC2CN1[*:1]"

  Scenario: Cartesian enumeration of the cell's core with R1 and R2 imported from a table gives 16 molecules
    When user picks "Enumerate Markush Structure..." from the context menu of the "cell 2 of Core" area of grid
    Then "Markush Enumerator" dialog should be visible
    And "Enumerator type" input in "Markush Enumerator" dialog should have the value "Cartesian"
    And OK button in "Markush Enumerator" dialog should be disabled
    And "Markush Enumerator" dialog should contain the text "R2"
    And "Markush Enumerator" dialog should not contain the text "R3"
    When user clicks on second "Import data" button in "Markush Enumerator" dialog
    Then "Import R-Groups" dialog should be visible
    When user selects "chem_enum_rgroups" in Table input in "Import R-Groups" dialog
    And user selects "R1" in Column input in "Import R-Groups" dialog
    Then "Target R#" input in "Import R-Groups" dialog should have the value "1"
    When user clicks on OK button in "Import R-Groups" dialog
    Then the "Import R-Groups" dialog should close
    And "Markush Enumerator" dialog should contain the text "r group 4 · R1"
    When user clicks on second "Import data" button in "Markush Enumerator" dialog
    Then "Import R-Groups" dialog should be visible
    When user selects "chem_enum_rgroups" in Table input in "Import R-Groups" dialog
    And user selects "R2" in Column input in "Import R-Groups" dialog
    And user enters "2" into "Target R#" input in "Import R-Groups" dialog
    And user clicks on OK button in "Import R-Groups" dialog
    Then the "Import R-Groups" dialog should close
    And "Markush Enumerator" dialog should contain the text "r group 4 · R2"
    And "Markush Enumerator" dialog should not contain the text "r group 5 · R1"
    And "Markush Enumerator" dialog should not contain the text "r group 5 · R2"
    And "Markush Enumerator" dialog should contain the text "16 molecules will be generated"
    And OK button in "Markush Enumerator" dialog should be enabled
    When user enters "BddMarkushCartesian" into "Table name" input in "Markush Enumerator" dialog
    And user clicks on OK button in "Markush Enumerator" dialog
    Then the "BddMarkushCartesian" view should be current
    And the table should have 16 rows
    And the table should have the columns "Enumerated, Core, R1, R2"
    And "Enumerated" column should have semantic type "Molecule"
    And "Core" column should have semantic type "Molecule"
    And "R1" column should have semantic type "Molecule"
    And "R2" column should have semantic type "Molecule"
    And "Enumerated" column should have no missing values
    And every value of "Core" column should match "^C1C2CN\(\[\*:2\]\)CC2CN1\[\*:1\]$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Zip enumeration of the same core and R-groups gives 4 molecules under the typed table name
    When user picks "Enumerate Markush Structure..." from the context menu of the "cell 2 of Core" area of grid
    Then "Markush Enumerator" dialog should be visible
    When user clicks on second "Import data" button in "Markush Enumerator" dialog
    And user selects "chem_enum_rgroups" in Table input in "Import R-Groups" dialog
    And user selects "R1" in Column input in "Import R-Groups" dialog
    And user clicks on OK button in "Import R-Groups" dialog
    Then the "Import R-Groups" dialog should close
    When user clicks on second "Import data" button in "Markush Enumerator" dialog
    And user selects "chem_enum_rgroups" in Table input in "Import R-Groups" dialog
    And user selects "R2" in Column input in "Import R-Groups" dialog
    And user enters "2" into "Target R#" input in "Import R-Groups" dialog
    And user clicks on OK button in "Import R-Groups" dialog
    Then the "Import R-Groups" dialog should close
    And "Markush Enumerator" dialog should contain the text "16 molecules will be generated"
    When user selects "Zip" in "Enumerator type" input in "Markush Enumerator" dialog
    Then "Markush Enumerator" dialog should contain the text "4 molecules will be generated"
    And "Markush Enumerator" dialog should not contain the text "16 molecules will be generated"
    When user enters "BddMarkushZip" into "Table name" input in "Markush Enumerator" dialog
    And user clicks on OK button in "Markush Enumerator" dialog
    Then the "BddMarkushZip" view should be current
    And the table should have 4 rows
    And every value of "Core" column should match "^C1C2CN\(\[\*:2\]\)CC2CN1\[\*:1\]$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Chem | Transform | Markush Enumeration opens the dialog on the current core
    When user clicks on the "cell 2 of Core" area of grid
    Then row 2 should be current
    When user picks "Chem > Transform > Markush Enumeration..." from the top menu
    Then "Markush Enumerator" dialog should be visible
    And "Markush Enumerator" dialog should contain the text "R2"
    And "Markush Enumerator" dialog should not contain the text "R3"
    When user clicks on CANCEL button in "Markush Enumerator" dialog
    Then the "Markush Enumerator" dialog should close
    And no errors should have been logged

  Scenario: The Markush Enumerator app keeps its content after the user leaves its view and comes back
    When user clicks on browse tab
    Given Apps tree node inside browse tree is expanded
    And Apps---Chem tree node inside browse tree is expanded
    When user double-clicks Apps---Chem---Markush-Enumerator tree node inside browse tree
    Then the "Markush Enumerator" view should be current
    And "Cores" text should be visible
    And "R-Groups" text should be visible
    And "Preview" text should be visible
    And ENUMERATE button should be visible
    And "Enumerator type" input should be visible
    When user clicks on ENUMERATE button
    Then "Enumerate" dialog should be visible
    And Output input in "Enumerate" dialog should have the value "New table"
    And "Table name" input in "Enumerate" dialog should be visible
    And "Remove duplicates" checkbox in "Enumerate" dialog should be visible
    And OK button in "Enumerate" dialog should be visible
    When user clicks on CANCEL button in "Enumerate" dialog
    Then the "Enumerate" dialog should close
    When user clicks on the tab of the "chem_enum_cores" view
    Then the "chem_enum_cores" view should be current
    And "Cores" text should be hidden
    When user clicks on the tab of the "Markush Enumerator" view
    Then the "Markush Enumerator" view should be current
    And "Cores" text should be visible
    And "R-Groups" text should be visible
    And ENUMERATE button should be visible
    And no errors should have been logged
