@journey @viewers @realizes:biostructureviewer.viewer.ngl
Feature: NGL viewer on a structure table — empty state and settings
  The NGL viewer does not take its structure from a Molecule3D column: added from the Add Viewer
  gallery to pdb_data.csv it shows its Open... link and no structure. Its Style > Representation
  offers the NGL representations. Translated from the TestTrack case
  BiostructureViewer/ngl-viewer-extension.

  Not translated: Scenario 1's entry, Show > NGL from a Molecule3D cell. That command docks the
  viewer's root straight into the view: the panel has no viewer name, the table view does not list
  the viewer, and neither its settings icon nor its canvas can be named by a step (requested in the
  request document, with the scenario). The settings claims of Scenario 1 are made here, on the
  viewer from the gallery: the settings do not depend on how the viewer was added. What the canvas
  draws in the chosen style is not a claim.

  Adding the viewer is checked for errors and balloons at the end of the Background, once the
  Open... link shows the viewer has come up empty.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree
    Then the "pdb_data" table view should open with 6 rows
    And "Add viewer" icon in toolbar should be visible
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user clicks on first "NGL" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And NGL viewer should be visible
    And "Open..." link in NGL viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: NGL added to a table with only a Molecule3D column shows Open... and no structure
    Then "Open..." link in NGL viewer should be visible
    And no errors should have been logged

  Scenario: Style > Representation offers the NGL representations and takes ball+stick
    When user clicks on settings icon of NGL viewer
    Given "Style" category in context panel is expanded
    Then "Representation" property in context panel should contain text "cartoon"
    When user clicks on "Representation" property in context panel
    Then "Representation" property in context panel should offer "cartoon, backbone, ball+stick, licorice, hyperball, surface"
    When user selects "ball+stick" in "Representation" property in context panel
    Then "representation" property of NGL viewer should be "ball+stick"
    When user selects "cartoon" in "Representation" property in context panel
    Then "representation" property of NGL viewer should be "cartoon"
    And no errors should have been logged
    And no error or warning balloon should have been shown
