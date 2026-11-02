Feature: A molecule column exported to CSV as SMILES
  Save > As CSV (options)... writes the selected columns only, and with Molecules as SMILES the
  MOLBLOCK column goes out as SMILES text: spgi's Structure is MOLBLOCK in all 100 rows, and the file
  holds no "M  END". Translated from TestTrack General/molecule-in-exported-csv.md. That the SMILES read
  back as a Molecule column, and that a plain export keeps the MOLBLOCKs, is the conversion's claim,
  tested in Chem (src/tests/save-as-csv-tests.ts) rather than through a download and an upload here.
  An empty molecule cell goes out as an empty field, with no placeholder text.

  The conversion calls Chem:convertNotation, so Chem must be installed. Nothing is put on the server;
  the dialog keeps its last options in the browser, so the scenario sets every box itself.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And user opens spgi dataset
    When user clicks on the "header Id" area of grid holding Control
    And user clicks on the "header Structure" area of grid holding Control
    And user clicks on the "header Chemist" area of grid holding Control
    Then columns "Id, Structure, Chemist" should be selected

  Scenario: Selected columns exported with Molecules as SMILES, and an empty molecule cell exported as an empty field
    Given user watches downloads
    When user clicks on Export icon in toolbar
    And user clicks on "As CSV (options)..." text in toolbar
    Then "Save as CSV" dialog should be visible
    When user checks "Molecules as Smiles" checkbox in "Save as CSV" dialog
    And user checks "Selected Columns Only" checkbox in "Save as CSV" dialog
    And user unchecks "Selected Rows Only" checkbox in "Save as CSV" dialog
    And user unchecks "Filtered Rows Only" checkbox in "Save as CSV" dialog
    And user downloads a file through OK button in "Save as CSV" dialog
    Then a file "spgi-100.csv" should have been downloaded
    And the downloaded file "spgi-100.csv" should contain text "Id,Structure,Chemist"
    And the downloaded file "spgi-100.csv" should contain 100 occurrences of "CAST-"
    And the downloaded file should not contain "M  END"
    And the downloaded file should not contain "Last Published Date"
    And the downloaded file should not contain "CAST Idea ID"
    And no errors should have been logged
    And no error or warning balloon should have been shown
    Given user opens a table "mols" with:
      | smiles   | id |
      | CCO      | 1  |
      |          | 2  |
      | c1ccccc1 | 3  |
    Then "smiles" column should have semantic type "Molecule"
    Given user watches downloads
    When user clicks on Export icon in toolbar
    And user clicks on "As CSV (options)..." text in toolbar
    Then "Save as CSV" dialog should be visible
    When user checks "Molecules as Smiles" checkbox in "Save as CSV" dialog
    And user unchecks "Selected Columns Only" checkbox in "Save as CSV" dialog
    And user unchecks "Selected Rows Only" checkbox in "Save as CSV" dialog
    And user unchecks "Filtered Rows Only" checkbox in "Save as CSV" dialog
    And user downloads a file through OK button in "Save as CSV" dialog
    Then a file "mols.csv" should have been downloaded
    And the downloaded file "mols.csv" should contain text ",2"
    And the downloaded file should not contain "null"
    And the downloaded file "mols.csv" should contain text "CCO,1"
    And no errors should have been logged
