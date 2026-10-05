@journey @viewers @realizes:biostructureviewer.viewer.biostructure
Feature: Biostructure viewer on a structure table — structure column, representation, camera, download
  The Biostructure viewer (the Mol* engine) added from the Add Viewer gallery to pdb_data.csv, the
  package's own table of six PDB structures (`pdb` holds full PDB texts, detected as Molecule3D;
  `pdb_id` holds PDB IDs). Nothing here reaches an outside service. Translated from the TestTrack
  case BiostructureViewer/biostructure-viewer (GROK-11759 for the representation switch).

  The md expects the viewer to take the `pdb` column by itself. On the stand it does not: the
  viewer added from the gallery keeps Biostructure Id empty and shows its Data file input
  (GROK-21119; its known-failure scenario is in file-open-and-preview). The Background therefore
  picks `pdb` in Biostructure Id through the settings, as a user facing the empty viewer does; what
  the viewer then shows is the subject of the scenarios.

  That the Mol* engine is built is claimed through its Reset Camera overlay button, which only a
  built engine has. What the 3D viewport draws (the representation, the camera) is not claimed:
  there is no reading of it, and pixels are not a claim. That a click on row 3 reaches the viewer
  is claimed by the download: after it, the viewer exports row 3's structure (2BDJ) — in PDB its
  CRYST1 record, 42.100 54.600 69.000 (row 1's reads 62.800 62.800 83.500), in CIF its polymer
  entity, PROTO-ONCOGENE TYROSINE-PROTEIN KINASE SRC, which no other row of the table holds. Mol*
  writes the PDB file without the source's HEADER records.

  Adding the viewer and building its engine are checked for errors and balloons at the end of the
  Background, after the Reset Camera button shows the engine is built. The representation a switch
  applies has no reading of its own (requested); each switch's error check comes after the
  property reads back and the engine's Reset Camera button is on the page.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree
    Then the "pdb_data" table view should open with 6 rows
    And "pdb" column should have semantic type "Molecule3D"
    And "Add viewer" icon in toolbar should be visible
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And Biostructure viewer should be visible
    When user clicks on settings icon of Biostructure viewer
    And user selects "pdb" in "Biostructure Id" property in context panel
    Then "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The viewer shows the pdb column's structure with the cartoon representation
    Then "Data File" input in Biostructure viewer should be absent
    Given "Style" category in context panel is expanded
    Then "Representation" property in context panel should contain text "cartoon"
    When user clicks on the "cell 3 of id" area of grid
    Then row 3 should be current
    And "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The representation switches from the settings (GROK-11759)
    When user selects "ball-and-stick" in "Representation" property in context panel
    Then "representation" property of Biostructure viewer should be "ball-and-stick"
    And "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user selects "molecular-surface" in "Representation" property in context panel
    Then "representation" property of Biostructure viewer should be "molecular-surface"
    And "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user selects "cartoon" in "Representation" property in context panel
    Then "representation" property of Biostructure viewer should be "cartoon"
    And "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Reset Camera overlay button keeps the viewer working
    When user clicks on "Reset Camera" button in Biostructure viewer
    Then "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The viewer's context menu downloads the structure as PDB and as CIF
    Given user watches downloads
    When user opens the context menu of Biostructure viewer
    Then the open menu should list "Download > As PDB"
    And the open menu should list "Download > As CIF"
    When user closes the context menu
    And user picks "Download > As PDB" from the context menu of Biostructure viewer
    Then a file "file.pdb" should have been downloaded
    And the downloaded file "file.pdb" should contain text "CRYST1   42.100   54.600   69.000"
    When user picks "Download > As CIF" from the context menu of Biostructure viewer
    Then a file "file.cif" should have been downloaded
    And the downloaded file "file.cif" should contain text "_atom_site"
    And the downloaded file "file.cif" should contain text "PROTO-ONCOGENE TYROSINE-PROTEIN KINASE SRC"
    And no errors should have been logged
    And no error or warning balloon should have been shown
