@journey @realizes:chem.cp.import-export-formats
Feature: Saving a table as SDF
  Opening an SDF or a MOL2 file is the importers' work, tested in Chem src/tests/detector-tests.ts
  (detectMolblockSDF) and src/tests/mol2-importer-tests.ts ('mol2 to SDF').

  The export icon of the toolbar offers As SDF..., which on smiles opens the Save as SDF dialog on canonical_smiles and downloads
  smiles.sdf, a V2000 molblock per row with the record terminator and the other columns as SD
  fields. Choosing the v3Kmolblock notation writes V3000 molblocks instead, and Filtered Rows Only
  writes one record for each of the 254 rows with two aromatic rings that pass the filter.

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: Save as SDF offers the molecule column and downloads a record per row
    And user opens smiles dataset
    And user watches downloads
    When user clicks on "arrow-to-bottom" icon in toolbar
    And user clicks on "As SDF..." item
    Then "Save as SDF" dialog should be visible
    And Molecules input in "Save as SDF" dialog should contain text "canonical_smiles"
    And Notation input in "Save as SDF" dialog should have value ""
    And "Visible Columns Only" input in "Save as SDF" dialog should be checked
    And "Selected Columns Only" input in "Save as SDF" dialog should not be checked
    And "Filtered Rows Only" input in "Save as SDF" dialog should not be checked
    And "Selected Rows Only" input in "Save as SDF" dialog should not be checked
    When user clicks on OK button in "Save as SDF" dialog
    Then the "Save as SDF" dialog should close
    And a file "smiles.sdf" should have been downloaded
    And the downloaded file "smiles.sdf" should contain text "M  END"
    And the downloaded file "smiles.sdf" should contain text "V2000"
    And the downloaded file "smiles.sdf" should contain text "$$$$"
    And the downloaded file "smiles.sdf" should contain text ">  <molregno>"
    And the downloaded file "smiles.sdf" should contain 1000 occurrences of "$$$$"
    And the table should have 1000 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Save as SDF writes V3000 molblocks when that notation is chosen
    And user opens smiles dataset
    And user watches downloads
    When user clicks on "arrow-to-bottom" icon in toolbar
    And user clicks on "As SDF..." item
    And user selects "v3Kmolblock" in Notation input in "Save as SDF" dialog
    And user clicks on OK button in "Save as SDF" dialog
    Then the "Save as SDF" dialog should close
    And a file "smiles.sdf" should have been downloaded
    And the downloaded file "smiles.sdf" should contain text "V3000"
    And the downloaded file "smiles.sdf" should contain text "M  V30 BEGIN ATOM"
    And the downloaded file "smiles.sdf" should contain 1000 occurrences of "$$$$"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Save as SDF writes only the filtered rows when asked to
    And user opens smiles dataset
    And user watches downloads
    When user filters rows where "NumAromaticRings" is between 2 and 2
    Then 254 rows should pass the filter
    When user clicks on "arrow-to-bottom" icon in toolbar
    And user clicks on "As SDF..." item
    And user checks "Filtered Rows Only" input in "Save as SDF" dialog
    And user clicks on OK button in "Save as SDF" dialog
    Then the "Save as SDF" dialog should close
    And a file "smiles.sdf" should have been downloaded
    And the downloaded file "smiles.sdf" should contain 254 occurrences of "$$$$"
    And no errors should have been logged
    And no error or warning balloon should have been shown
