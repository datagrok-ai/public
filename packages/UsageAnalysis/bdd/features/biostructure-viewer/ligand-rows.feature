@journey @viewers @realizes:biostructureviewer.viewer.biostructure @realizes:biostructureviewer.viewer.ngl
Feature: Which rows' ligands the Biostructure and NGL viewers show — selected, current and mouse-over rows
  Both viewers on 1bdq-obs-pred.sdf (22 ligand poses) over the protein 1bdq.pdb, from App Data >
  BiostructureViewer > samples: the Biostructure viewer gets the protein through its Data file input,
  the NGL viewer through its Open... link. The Behaviour category decides whose ligands a viewer
  draws: the selected rows (Show Selected Rows Ligands, off by default), the current row and the row
  under the pointer (both on by default). Translated from the TestTrack cases
  BiostructureViewer/ngl-viewer-extension (steps 7-9) and BiostructureViewer/property-surface-extension
  (scenario 4), and from the GROK-17967 parity spec: the selected rows' ligands must be the same in
  both engines.

  Which rows' ligands a viewer drew is its `ligand rows` reading (distinct 1-based rows) and how many
  its `ligands shown`. The current and mouse-over sources are switched off before the selection is
  made, so the four ligands can only come from the selection. A source that is off is claimed by the
  reading staying as it was once the viewer has nothing pending. That the default current-row source
  works is in data-file-persistence.feature; here the NGL viewer's two toggles are claimed both ways.
  What the 3D canvas draws is not a claim.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree
    Then the "1bdq-obs-pred" table view should open with 22 rows
    When user clicks on "Add viewer" icon in toolbar
    And user types "Biostructure" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then "Data File" input in Biostructure viewer should be visible
    When user clicks on "folder-tree" icon in "Data File" input in Biostructure viewer
    Then "Select a file" dialog should be visible
    Given Files---App-Data tree node in "Select a file" dialog is expanded
    And Files---App-Data---BiostructureViewer tree node in "Select a file" dialog is expanded
    And Files---App-Data---BiostructureViewer---samples tree node in "Select a file" dialog is expanded
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node in "Select a file" dialog
    And user clicks on OK button in "Select a file" dialog
    Then "Reset Camera" button in Biostructure viewer should be visible
    And "ligandColumnName" property of Biostructure viewer should be "molecule"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Biostructure shows the ligands of the selected rows with Show Selected Rows Ligands
    When user clicks on settings icon of Biostructure viewer
    Given "Behaviour" category in context panel is expanded
    When user unchecks "Show Current Row Ligand" property in context panel
    And user unchecks "Show Mouse-Over Row Ligand" property in context panel
    Then the "ligands shown" reading of Biostructure viewer should be 0
    When user takes a snapshot of Biostructure viewer
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then rows 1 to 4 should be selected
    And 4 rows should be selected
    And the "ligands shown" reading of Biostructure viewer should be the same as before
    When user clicks on settings icon of Biostructure viewer
    And user checks "Show Selected Rows Ligands" property in context panel
    Then the "ligand rows" reading of Biostructure viewer should be "1, 2, 3, 4"
    And the "ligands shown" reading of Biostructure viewer should be 4
    When user presses Escape in grid
    Then no rows should be selected
    And the "ligands shown" reading of Biostructure viewer should be 0
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: NGL shows the same selected rows' ligands as Biostructure (GROK-17967)
    When user clicks on "Add viewer" icon in toolbar
    And user types "NGL" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "NGL" button in "Add Viewer" dialog
    Then NGL viewer should be visible
    When user uploads "../../BiostructureViewer/files/samples/1bdq.pdb" through "Open..." link in NGL viewer
    Then the "structure loaded" reading of NGL viewer should be "true"
    When user clicks on settings icon of NGL viewer
    Given "Data" category in context panel is expanded
    When user selects "molecule" in "Ligand" property in context panel
    Then "ligandColumnName" property of NGL viewer should be "molecule"
    Given "Behaviour" category in context panel is expanded
    When user unchecks "Show Current Row Ligand" property in context panel
    And user unchecks "Show Mouse-Over Row Ligand" property in context panel
    And user checks "Show Selected Rows Ligands" property in context panel
    Then the "ligands shown" reading of NGL viewer should be 0
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 4 rows should be selected
    And the "ligand rows" reading of NGL viewer should be "1, 2, 3, 4"
    And the "ligands shown" reading of NGL viewer should be 4
    And the "ligand rows" reading of Biostructure viewer should be "1, 2, 3, 4"
    When user presses Escape in grid
    Then no rows should be selected
    And the "ligands shown" reading of NGL viewer should be 0
    And the "ligands shown" reading of Biostructure viewer should be 0
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: NGL's Show Current Row Ligand and Show Mouse-Over Row Ligand switch their ligand on and off
    When user clicks on settings icon of NGL viewer
    Given "Behaviour" category in context panel is expanded
    When user checks "Show Current Row Ligand" property in context panel
    And user clicks on the "cell 6 of molecule" area of grid
    Then row 6 should be current
    And the "ligand rows" reading of NGL viewer should be "6"
    And the "ligands shown" reading of NGL viewer should be 1
    When user clicks on settings icon of NGL viewer
    And user unchecks "Show Current Row Ligand" property in context panel
    Then the "ligands shown" reading of NGL viewer should be 0
    When user takes a snapshot of NGL viewer
    And user hovers over the "cell 9 of molecule" area of grid
    Then the mouse-over row of the table should be 9
    And the "ligands shown" reading of NGL viewer should be the same as before
    When user clicks on settings icon of NGL viewer
    And user checks "Show Mouse-Over Row Ligand" property in context panel
    And user hovers over the "cell 9 of molecule" area of grid
    Then the mouse-over row of the table should be 9
    And the "ligand rows" reading of NGL viewer should be "9"
    And the "ligands shown" reading of NGL viewer should be 1
    And no errors should have been logged
    And no error or warning balloon should have been shown
