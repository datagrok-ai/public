Feature: Merging a table from the user's files into a Hit Design campaign
  The ribbon's merge icon opens "Merge with table". Its OK stays disabled while the dialog lacks what the merge needs,
  and a line under the inputs says why: no source file yet, the V-iD column missing for a merge by V-iD (the
  default), the molecule column missing for a merge by canonical SMILES. By SMILES, a campaign row whose molecule is the file's (canonically: OCC is CCO) takes
  the file's values, and a molecule the campaign lacks comes in as a new row with the first stage, a V-iD registered
  for it and the template's compute (MW); by V-iD, the row holding the file's V-iD takes its values and nothing is
  added without a molecule column. The merged campaign is saved, with the balloon "Table merged into campaign.".

  The source files are CSVs in the user's home folder, typed into Source file as a share path is pasted or picked in
  the platform's "Select a file" dialog ("Add file from Files"); they go at feature end. A file picked from a share
  is remembered with the campaign (with Auto-merge on, the campaign reads and merges it again whenever it opens: the
  package test "Compat: campaign.json" pins that); a file uploaded from the computer is merged once. The campaign is
  made for the feature through the JS API (bindings/campaign.ts); its V-iDs are plain values for the merge to match,
  not dictionary entries. The new molecule is one only bdd features use, so the
  dictionary row its registration makes is the feature's own and goes with the campaign.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Merge" with the stages "Idea, Design" and a table of:
      | Molecule | Stage  | V-iD    | MW    |
      | CCO      | Design | V900001 | 46.07 |
      | c1ccccc1 | Design | V900002 | 78.11 |
    And the user's file "bdd-merge-by-smiles.csv" holds the table:
      | Molecule                   | Score |
      | OCC                        | 5     |
      | NCCCCCCCCCCCCCCCCCCCCCCCCO | 7     |
    And the user's file "bdd-merge-by-vid.csv" holds the table:
      | V-iD    | Priority |
      | V900002 | High     |
    When user opens the campaign "BDD Merge" by its link
    And user clicks on "Merge another table into the campaign" icon in toolbar
    Then "Merge with table" dialog should be visible

  Scenario: OK stays disabled, saying why, until the dialog has what the merge needs
    Then OK button in "Merge with table" dialog should be disabled
    And "Merge with table" dialog should contain text "Pick a source file to merge from."
    And CANCEL button in "Merge with table" dialog should be enabled
    When user enters the path of the user's file "bdd-merge-by-smiles.csv" into "Source file" input in "Merge with table" dialog
    Then "Merge with table" dialog should contain text "Pick the V-iD column from the incoming table."
    And OK button in "Merge with table" dialog should be disabled
    When user selects "By canonical SMILES" in "Merge by" input in "Merge with table" dialog
    Then OK button in "Merge with table" dialog should be enabled
    When user selects "(none)" in "Molecule column (incoming)" input in "Merge with table" dialog
    Then "Merge with table" dialog should contain text "Pick the molecule column from the incoming table."
    And OK button in "Merge with table" dialog should be disabled
    When user clicks on CANCEL button in "Merge with table" dialog
    Then "Merge with table" dialog should be absent
    And the table should have 2 rows

  Scenario: By canonical SMILES, the matching row takes the file's values and the new molecule comes in with a stage, a V-iD and MW
    When user enters the path of the user's file "bdd-merge-by-smiles.csv" into "Source file" input in "Merge with table" dialog
    And user selects "By canonical SMILES" in "Merge by" input in "Merge with table" dialog
    And user clicks on OK button in "Merge with table" dialog
    Then an info balloon containing "Table merged into campaign" should have been shown
    And the campaign "BDD Merge" should be saved with 3 rows
    And the table should have 3 rows
    And the value of "Score" column in row 1 should be "5"
    And the value of "Score" column in row 2 should be ""
    And the value of "Molecule" column in row 3 should be "NCCCCCCCCCCCCCCCCCCCCCCCCO"
    And the value of "Stage" column in row 3 should be "Idea"
    And the value of "V-iD" column in row 3 should match "^V[0-9]{6}$"
    And the value of "MW" column in row 3 should match "^369[.]6"
    And the campaign "BDD Merge" should be saved with "5" in row 1 of the "Score" column

  Scenario: By V-iD, the row of that V-iD takes the file's values and no row is added
    When user enters the path of the user's file "bdd-merge-by-vid.csv" into "Source file" input in "Merge with table" dialog
    Then "V-iD column (incoming)" input in "Merge with table" dialog should have value "V-iD"
    And OK button in "Merge with table" dialog should be enabled
    When user clicks on OK button in "Merge with table" dialog
    Then an info balloon containing "Table merged into campaign" should have been shown
    And the campaign "BDD Merge" should be saved with "High" in row 2 of the "Priority" column
    And the table should have 2 rows
    And the value of "Priority" column in row 2 should be "High"
    And the value of "Priority" column in row 1 should be ""

  # a file from Files: remembered with the campaign, to be merged on every open (package test "Compat:
  # campaign.json" pins the merge on open)
  Scenario: A file picked from Files is merged and remembered with the campaign, and the dialog opens on it again
    When user clicks on "Add file from Files" icon in "Source file" input in "Merge with table" dialog
    Then "Select a file" dialog should be visible
    Given "My files" tree node in "Select a file" dialog is expanded
    When user double-clicks on "bdd-merge-by-vid.csv" tree node in "Select a file" dialog
    Then "Select a file" dialog should be absent
    And "Source file" input in "Merge with table" dialog should hold the path of the user's file "bdd-merge-by-vid.csv"
    And "V-iD column (incoming)" input in "Merge with table" dialog should have value "V-iD"
    When user checks "Auto-merge on campaign open" input in "Merge with table" dialog
    And user clicks on OK button in "Merge with table" dialog
    Then an info balloon containing "Table merged into campaign" should have been shown
    And the campaign "BDD Merge" should be saved with "High" in row 2 of the "Priority" column
    And the campaign "BDD Merge" should be saved to merge the user's file "bdd-merge-by-vid.csv" on every open
    When user clicks on "Merge another table into the campaign" icon in toolbar
    Then "Source file" input in "Merge with table" dialog should hold the path of the user's file "bdd-merge-by-vid.csv"
    And "V-iD column (incoming)" input in "Merge with table" dialog should have value "V-iD"
    And "Auto-merge on campaign open" input in "Merge with table" dialog should be checked

  # a file from the computer: merged once, nothing to read again on open
  Scenario: A file uploaded from the computer is merged once, with Auto-merge on campaign open off
    When user clicks on "Open file" icon in "Source file" input in "Merge with table" dialog and chooses the "fixtures/merge-local.csv" file of the project
    Then "V-iD column (incoming)" input in "Merge with table" dialog should have value "V-iD"
    And "Auto-merge on campaign open" input in "Merge with table" dialog should be disabled
    When user clicks on OK button in "Merge with table" dialog
    Then an info balloon containing "Table merged into campaign" should have been shown
    And the campaign "BDD Merge" should be saved with "Urgent" in row 1 of the "Priority" column
    And the campaign "BDD Merge" should be saved to merge no file on open
    And the value of "Priority" column in row 1 should be "Urgent"
