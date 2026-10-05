Feature: A molecule column exported to CSV as SMILES reads back as molecules
  Save > As CSV (options)... writes the selected columns only, and with Molecules as SMILES the
  MOLBLOCK column goes out as SMILES text; the file opened again has the same columns and a
  Molecule column. spgi's Structure is MOLBLOCK in all 100 rows, so a plain export of the same
  selection carries "M  END" — the second scenario, which shows it is the option that converts.
  Translated from TestTrack General/molecule-in-exported-csv.md. The old spec never opened the dialog
  (it called toCsvEx and DG.Utils.download), so nothing of it is kept.

  The conversion calls Chem:convertNotation, so Chem must be installed. Nothing is put on the server;
  the dialog keeps its last options in the browser, so every scenario sets both boxes itself.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And user opens spgi dataset
    When user clicks on the "header Id" area of grid holding Control
    And user clicks on the "header Structure" area of grid holding Control
    And user clicks on the "header Chemist" area of grid holding Control
    Then columns "Id, Structure, Chemist" should be selected

  Scenario: Selected columns exported with Molecules as SMILES, and the file opened again
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
    When user closes all views
    And user uploads the downloaded file through "Open local file" icon inside browse toolbar
    Then the table should have 100 rows
    And the table should have 3 columns
    And the table should have the columns "Id, Structure, Chemist"
    And the value of "Id" column in row 1 should be "CAST-634783"
    And "Structure" column should have semantic type "Molecule"
    And "Structure" column should have units "smiles"
    And every value of "Structure" column should match "^\S+$"
    And "Structure" column should have no missing values
    And "Structure" column should have at least 100 distinct values
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Without the option the same selection carries the MOLBLOCKs
    Given user watches downloads
    When user clicks on Export icon in toolbar
    And user clicks on "As CSV (options)..." text in toolbar
    Then "Save as CSV" dialog should be visible
    When user unchecks "Molecules as Smiles" checkbox in "Save as CSV" dialog
    And user checks "Selected Columns Only" checkbox in "Save as CSV" dialog
    And user unchecks "Selected Rows Only" checkbox in "Save as CSV" dialog
    And user unchecks "Filtered Rows Only" checkbox in "Save as CSV" dialog
    And user downloads a file through OK button in "Save as CSV" dialog
    Then the downloaded file "spgi-100.csv" should contain text "Id,Structure,Chemist"
    And the downloaded file should contain "M  END"
    And the downloaded file should not contain "Last Published Date"
    And the downloaded file should not contain "CAST Idea ID"
    And no errors should have been logged
    And no error or warning balloon should have been shown
