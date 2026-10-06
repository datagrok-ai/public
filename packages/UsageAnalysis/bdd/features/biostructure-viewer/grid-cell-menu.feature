@journey @viewers @realizes:biostructureviewer.cell.molecule3d @realizes:biostructureviewer.viewer.biostructure
Feature: Grid context menu on Molecule3D cells (GROK-14552)
  On a Molecule3D cell the package adds Copy, Download and a Show group with Biostructure and NGL
  to the grid's context menu; a PDB_ID cell gets none of them. Translated from the TestTrack case
  BiostructureViewer/biostructureviewer-bug-grok-14552, on pdb_data.csv (`pdb` is Molecule3D,
  `pdb_id` is PDB_ID).

  Show > Biostructure and Show > NGL dock a viewer of the table view, which is claimed visible with
  the structure it loaded (`structure loaded`) and closed again. The right-click on the empty area
  of a row past the last column, the bug of GROK-14552 itself, comes last (pdb_data's columns end
  well short of the grid's right edge).

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
    And "pdb_id" column should have semantic type "PDB_ID"

  Scenario: A Molecule3D cell's menu offers Copy, Download and Show > Biostructure / NGL
    When user right-clicks on the "cell 1 of pdb" area of grid
    Then the open menu should list "Copy"
    And the open menu should list "Download"
    And the open menu should list "Show > Biostructure"
    And the open menu should list "Show > NGL"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Copy puts the cell's PDB text on the clipboard
    When user picks "Copy" from the context menu of the "cell 1 of pdb" area of grid
    Then an info balloon containing "Value copied to clipboard" should have been shown
    And the clipboard should contain text "HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY"
    And no errors should have been logged

  Scenario: Show > Biostructure docks a Biostructure viewer with the cell's structure
    Then Biostructure viewer should be absent
    When user picks "Show > Biostructure" from the context menu of the "cell 2 of pdb" area of grid
    Then Biostructure viewer should be visible
    And the "structure loaded" reading of Biostructure viewer should be "true"
    And "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on close icon of Biostructure viewer
    Then Biostructure viewer should be absent

  Scenario: Show > NGL docks an NGL viewer with the cell's structure
    Then NGL viewer should be absent
    When user picks "Show > NGL" from the context menu of the "cell 2 of pdb" area of grid
    Then NGL viewer should be visible
    And the "structure loaded" reading of NGL viewer should be "true"
    And "Open..." link in NGL viewer should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on close icon of NGL viewer
    Then NGL viewer should be absent

  Scenario: A PDB_ID cell's menu has none of the Molecule3D items
    When user right-clicks on the "cell 1 of pdb_id" area of grid
    Then the open menu should list "Properties..."
    And the open menu should not list "Show"
    And the open menu should not list "Copy"
    And the open menu should not list "Download"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Right-clicking the row's empty space past the last column raises nothing (GROK-14552)
    When user right-clicks on empty space of row 2 of grid
    Then the open menu should list "Properties..."
    And the open menu should not list "Show"
    When user closes the context menu
    Then no errors should have been logged
    And no error or warning balloon should have been shown
