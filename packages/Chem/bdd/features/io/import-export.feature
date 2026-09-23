@journey @realizes:chem.cp.import-export-formats
Feature: Opening chemical file formats and saving a table as SDF
  Each chemical format opens through its own Chem importer. mol1K.sdf gives 1000 rows in a molecule
  column of V2000 molblocks, typed Molecule, with the five SD fields beside it; molecules.mol2 gives
  three rows in a molecules column, one per TRIPOS block. The cells of
  every one of them are drawn as coloured structures, while a number column of the same grid carries
  no colour at all.

  The export icon of the toolbar offers As SDF..., which on smiles opens the Save as SDF dialog on canonical_smiles and downloads
  smiles.sdf, a V2000 molblock per row with the record terminator and the other columns as SD
  fields. Choosing the v3Kmolblock notation writes V3000 molblocks instead.

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: An SDF file opens as one row per record, with its fields beside the molecule
    Given user opens mol1K.sdf dataset
    Then the table should have 1000 rows
    And the table should have a column "molecule"
    And "molecule" column should have semantic type "Molecule"
    And "molecule" column should have units "molblock"
    And "molecule" column should have no missing values
    And every value of "molecule" column should contain "V2000"
    And every value of "molecule" column should contain "M  END"
    And the table should have a column "prID"
    And the table should have a column "pIC50_HIV_Integrase"
    And the table should have a column "Activity_Integrase"
    And the "cell 1 of molecule" area of grid should be painted
    And the "cell 1 of molecule" area of grid should be painted in at least 2 colors
    And the "cell 2 of molecule" area of grid should be painted in at least 2 colors
    And the "cell 1 of prID" area of grid should not contain the color "#FF0000"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A MOL2 file opens as one row per TRIPOS molecule
    Given user opens molecules.mol2 dataset
    Then the table should have 3 rows
    And the table should have 1 column
    And the table should have a column "molecules"
    And "molecules" column should have semantic type "Molecule"
    And "molecules" column should have units "molblock"
    And "molecules" column should have no missing values
    And every value of "molecules" column should contain "M  END"
    And the "cell 1 of molecules" area of grid should be painted
    And the "cell 1 of molecules" area of grid should be painted in at least 2 colors
    And no errors should have been logged
    And no error or warning balloon should have been shown

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
    Then fewer than 1000 rows should pass the filter
    When user clicks on "arrow-to-bottom" icon in toolbar
    And user clicks on "As SDF..." item
    And user checks "Filtered Rows Only" input in "Save as SDF" dialog
    And user clicks on OK button in "Save as SDF" dialog
    Then the "Save as SDF" dialog should close
    And a file "smiles.sdf" should have been downloaded
    And the downloaded file "smiles.sdf" should contain fewer than 1000 occurrences of "$$$$"
    And no errors should have been logged
    And no error or warning balloon should have been shown
