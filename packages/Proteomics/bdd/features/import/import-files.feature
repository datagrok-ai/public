Feature: Proteomics > Import opens a MaxQuant proteinGroups.txt as a table of protein groups, directly or through the column mapping
  The Proteomics menu's Import group reads a search engine's output from a file the user chooses
  (package README, "Supported inputs" and "Import": MaxQuant `proteinGroups.txt`, "contaminant &
  reverse-hit filtering"; the import is "also available from the left-bar main menu with no table
  open"). MaxQuant... opens the file as a new table view named after it, with the contaminant and
  reverse hits left out, the protein IDs typed as protein identifiers, and an info balloon that counts the
  imported protein groups. Generic Matrix... reads any CSV/TSV into its column-mapping dialog ("Import Generic
  Matrix": "Column-mapping dialog; auto-detects delimiter and log2 status"), which proposes the protein ID,
  gene name and intensity columns and previews them; OK opens the file as is, every row kept.

  The file is the package's own `files/demo/proteinGroups.txt` (27 protein groups, of which P04264 and
  CON__P02768 are marked "Potential contaminant" and REV__Q9Y6K1 "Reverse"), chosen in the browser's file
  chooser the menu opens. Left out on purpose:
  - the other file importers (Spectronaut Report and Candidates, FragPipe): the same file chooser path,
    their parsing is the package's parser tests (`src/tests`);
  - the log2 transform of the intensities and the gene-label resolution: numbers the parser tests check.
  Nothing is put on the server.

  Scenario: MaxQuant... opens the chosen proteinGroups.txt without its contaminant and reverse hits
    Given user is logged in
    When user clicks on Menu tab
    And user picks "Proteomics > Import > MaxQuant..." from the open menu and chooses the "../files/demo/proteinGroups.txt" file of the project
    Then the "proteinGroups" table view should open with 24 rows
    And an info balloon containing "Imported 24 protein groups" should have been shown
    And the value of "Protein IDs" column in row 1 should be "P04406"
    And "Protein IDs" column should have semantic type "Proteomics-ProteinId"
    And every value of "Protein IDs" column should match "^[A-Z][0-9A-Z]{5}$"
    And no errors should have been logged

  Scenario: Generic Matrix... maps the chosen file's columns in a dialog, and OK opens every row
    Given user is logged in
    When user clicks on Menu tab
    And user picks "Proteomics > Import > Generic Matrix..." from the open menu and chooses the "../files/demo/proteinGroups.txt" file of the project
    Then "Import Generic Matrix" dialog should be visible
    And editor of "Protein ID" input in "Import Generic Matrix" dialog should have text "Protein IDs"
    And editor of "Gene Name (optional)" input in "Import Generic Matrix" dialog should have text "Gene names"
    And "Intensity Columns" input in "Import Generic Matrix" dialog should contain text "(6) LFQ intensity Sample1"
    And "Log2 Transform" checkbox in "Import Generic Matrix" dialog should be checked
    And "Import Generic Matrix" dialog should contain text "Data appears to be raw intensities"
    And "Import Generic Matrix" dialog should contain text "Preview: 8 columns, 5 rows"
    When user clicks on OK button in "Import Generic Matrix" dialog
    Then the "proteinGroups" table view should open with 27 rows
    And an info balloon containing "Imported 27 proteins from proteinGroups.txt" should have been shown
    And "Protein IDs" column should have semantic type "Proteomics-ProteinId"
    And no errors should have been logged
