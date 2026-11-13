@guide @demo
Feature: Bioinformatics · Protein structures
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). BiostructureViewer previews structure files with the Mol* engine
  straight from Browse > Files > App Data > BiostructureViewer > samples (a PDB, then an mmCIF), and opens
  pdb_data.csv, six PDB entries whose pdb column holds whole PDB texts (Molecule3D); the Biostructure viewer,
  added from the Add Viewer gallery, shows the structure of the current row and follows it. Off camera the
  grid keeps the PDB id, the class and the description, in rows under 40 px: into a taller cell the PDB id
  renderer draws RCSB's picture of the entry, fetched from cdn.rcsb.org, and the clip shows nothing an
  outside service sends. Nothing is made on the server.

  Scenario: Protein structures in 3D
    Structure files preview in Mol* straight from Browse, and a 3D viewer follows a table of structures row by row.

    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree
    Then the "1bdq.pdb" view should be current
    And "Reset Camera" button should be visible
    When user clicks on Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree
    Then the "1RQ9.mmcif" view should be current
    And "Reset Camera" button should be visible
    When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree
    Then the "pdb_data" table view should open with 6 rows
    And "pdb" column should have semantic type "Molecule3D"
    Given the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "pdb_id, prop2, desc"
    And the rows of the grid are 32 px tall
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "Biostructure" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And Biostructure viewer should be visible
    Given the viewers have settled
    When user clicks on the "cell 3 of pdb_id" area of grid
    Then row 3 should be current
    Given the viewers have settled
    When the video says "The viewer follows the current row"
    And user clicks on the "cell 6 of pdb_id" area of grid
    Then row 6 should be current
    Given the viewers have settled
    Then no errors should have been logged
