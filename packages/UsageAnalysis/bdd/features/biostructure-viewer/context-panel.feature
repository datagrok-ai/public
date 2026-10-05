@journey @viewers @realizes:biostructureviewer.panel.3d-structure @realizes:biostructureviewer.panel.pdb-information
Feature: Context panel of a Molecule3D cell — 3D Structure and PDB Information
  The package's panes for the current Molecule3D cell: 3D Structure, a Biostructure viewer embedded
  in the pane, and PDB Information, header data read from the PDB text itself (no network call).
  Translated from the TestTrack case BiostructureViewer/context-panel-widgets-extension, on
  pdb_data.csv (row 1 is 1QBS, row 2 is 1ZP8).

  The embedded viewer is claimed by the Reset Camera button its Mol* engine builds inside the pane.
  The loader shown while the pane is rebuilt for the second cell is not readable. That the panel
  was rebuilt for row 2 is shown by PDB Information, whose values differ between the rows: the 3D
  Structure claim for row 2 is made only after PDB Information reads row 2's header.

  The Protein-Ligand Interactions and PDB id viewer panes (a server-side Python script, RCSB) are
  never expanded here: a pane's expanded state persists for the page, and both are manual-only
  (biostructureviewer-network-ui). PDB Information stays expanded from the first scenario into the
  second: a click on a pane header made right after the click on the next cell keeps the context
  panel on the previous cell (a click faster than a user makes), so the second scenario
  only reads the pane after its row changes. The last scenario collapses what is still expanded.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree
    Then the "pdb_data" table view should open with 6 rows
    Given the context panel is open

  Scenario: 3D Structure embeds a Biostructure viewer for the current Molecule3D cell
    When user clicks on the "cell 1 of pdb" area of grid
    Then row 1 should be current
    Given "3D Structure" section in context panel is expanded
    Then "Reset Camera" button in "3D Structure" section in context panel should be visible
    Given "PDB Information" section in context panel is expanded
    Then "PDB Information" section in context panel should contain text "ASPARTYL PROTEASE"
    When user clicks on the "cell 2 of pdb" area of grid
    Then row 2 should be current
    And "PDB Information" section in context panel should contain text "HYDROLASE"
    And "PDB Information" section in context panel should not contain text "ASPARTYL PROTEASE"
    And "Reset Camera" button in "3D Structure" section in context panel should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user collapses "3D Structure" section in context panel
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: PDB Information reads the header of the current cell's PDB text
    When user clicks on the "cell 1 of pdb" area of grid
    Then row 1 should be current
    And "PDB Information" section in context panel should contain text "ASPARTYL PROTEASE"
    And "PDB Information" section in context panel should contain text "HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY"
    And "PDB Information" section in context panel should contain text "rcsb.org/structure/1QBS"
    When user clicks on the "cell 2 of pdb" area of grid
    Then row 2 should be current
    And "PDB Information" section in context panel should contain text "HYDROLASE"
    And "PDB Information" section in context panel should contain text "HIV PROTEASE WITH INHIBITOR AB-2"
    And "PDB Information" section in context panel should contain text "rcsb.org/structure/1ZP8"
    And "PDB Information" section in context panel should not contain text "ASPARTYL PROTEASE"
    When user collapses "PDB Information" section in context panel
    Then no errors should have been logged
    And no error or warning balloon should have been shown
