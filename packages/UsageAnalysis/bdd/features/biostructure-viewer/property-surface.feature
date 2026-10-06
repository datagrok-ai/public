@journey @viewers @realizes:biostructureviewer.viewer.biostructure
Feature: Biostructure viewer settings — Binding Site Whole Residues and the Controls category
  Two parts of the viewer's settings no other feature exercises: the Binding Site Whole Residues
  switch and the Controls category. Translated from the TestTrack case
  BiostructureViewer/property-surface-extension.

  Whole residues against single atoms differs only in the drawing, which is not a claim; what is
  claimed is that the switch reaches the viewer and that the Binding Site settings reach the
  overlay's Binding site popover (its "Show side chains" box follows Show Binding Site).

  What Mol* shows when Show Import Controls is on is not checked (it is up to Mol*). Only the
  category's two switches and their defaults are claimed.

  The Background picks the `pdb` column in Biostructure Id through the settings, because the viewer
  does not take it by itself on the stand (GROK-21119). Adding the
  viewer and building its engine are checked for errors and balloons at the end of the Background.

  Before the viewer is added, a click on a pdb cell puts the table's own cell into the context
  panel. Run after the NGL feature, whose last object in the panel is an NGL viewer that is closed
  with its view, the settings click on the new viewer otherwise left the panel on the closed NGL
  viewer's settings in 8 runs out of 14 (the click comes faster than a user makes it).

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
    When user clicks on the "cell 1 of pdb" area of grid
    Then "PDB Information" section in context panel should be visible
    And "Add viewer" icon in toolbar should be visible
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "Biostructure" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And Biostructure viewer should be visible
    When user clicks on settings icon of Biostructure viewer
    And user selects "pdb" in "Biostructure Id" property in context panel
    Then "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Binding Site Whole Residues switches off and on under Show Binding Site
    Given "Binding Site" category in context panel is expanded
    Then "Show Binding Site" property in context panel should be unchecked
    And "Binding Site Whole Residues" property in context panel should be checked
    When user checks "Show Binding Site" property in context panel
    Then "showBindingSite" property of Biostructure viewer should be "true"
    And "Binding site" overlay button in Biostructure viewer should be selected
    When user clicks on "Binding site" overlay button in Biostructure viewer
    Then "Show side chains" checkbox should be visible
    And "Show side chains" checkbox should be checked
    When user presses Escape
    Then "Show side chains" checkbox should be hidden
    When user unchecks "Binding Site Whole Residues" property in context panel
    Then "bindingSiteWholeResidues" property of Biostructure viewer should be "false"
    And "Reset Camera" button in Biostructure viewer should be visible
    When user checks "Binding Site Whole Residues" property in context panel
    Then "bindingSiteWholeResidues" property of Biostructure viewer should be "true"
    When user unchecks "Show Binding Site" property in context panel
    Then "showBindingSite" property of Biostructure viewer should be "false"
    When user clicks on "Binding site" overlay button in Biostructure viewer
    Then "Show side chains" checkbox should be visible
    And "Show side chains" checkbox should be unchecked
    When user presses Escape
    Then "Show side chains" checkbox should be hidden
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Controls category holds Show Welcome Toast and Show Import Controls, both off
    Given "Controls" category in context panel is expanded
    Then "Show Welcome Toast" property in context panel should be unchecked
    And "Show Import Controls" property in context panel should be unchecked
    And no errors should have been logged
