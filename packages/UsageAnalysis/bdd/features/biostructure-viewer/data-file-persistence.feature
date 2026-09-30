@journey @serial @viewers @realizes:biostructureviewer.viewer.biostructure @realizes:views.projects
Feature: A structure loaded into an empty Biostructure viewer survives saving and reopening the project
  A Biostructure viewer added to a table with no structure column shows a Data file input
  (GROK-16143); the file chosen there is kept in the viewer's dataJson, so it must come back when
  the project is saved and reopened (GROK-17485). Translated from the TestTrack case
  BiostructureViewer/biostructureviewer-data-file-persistence, on 1bdq-obs-pred.sdf (22 ligand
  poses) and the protein 1bdq.pdb, both in App Data > BiostructureViewer > samples.

  The Mol* engine is claimed by its Reset Camera overlay button; the empty state by the Data file
  input. Kept without one line of Scenario 1: that Ligand Column Name names the molecule column —
  on the stand the viewer leaves it empty (the same suspected defect as the structure column not
  being taken, in the request document).

  Not translated: Scenario 2 (GROK-17967, only the current row's ligand is shown, in Mol* and in
  NGL): how many structures a viewer has loaded is not on the page — a package reading is
  requested, with the scenario, in the request document.

  The project has a fixed name and is removed with its tables when the feature starts and ends;
  serial, as every feature that saves a project and searches the Dashboards gallery.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And no project named "BsvDataFilePersistence" is on the server
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded

  Scenario: The empty Biostructure viewer loads a structure chosen in its Data file input (GROK-16143)
    When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree
    Then the "1bdq-obs-pred" table view should open with 22 rows
    And "molecule" column should have semantic type "Molecule"
    And "Add viewer" icon in toolbar should be visible
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And "Data File" input in Biostructure viewer should be visible
    And "Reset Camera" button in Biostructure viewer should be absent
    When user clicks on "folder-tree" icon in "Data File" input in Biostructure viewer
    Then "Select a file" dialog should be visible
    Given Files---App-Data tree node in "Select a file" dialog is expanded
    And Files---App-Data---BiostructureViewer tree node in "Select a file" dialog is expanded
    And Files---App-Data---BiostructureViewer---samples tree node in "Select a file" dialog is expanded
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node in "Select a file" dialog
    And user clicks on OK button in "Select a file" dialog
    Then "Reset Camera" button in Biostructure viewer should be visible
    And "Data File" input in Biostructure viewer should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The structure comes back when the project is saved and reopened (GROK-17485)
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BsvDataFilePersistence" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BsvDataFilePersistence" uploaded' should have been shown
    And "Share BsvDataFilePersistence" dialog should be visible
    When user clicks on CANCEL button in "Share BsvDataFilePersistence" dialog
    Then the "Share BsvDataFilePersistence" dialog should close
    And 1 project named "BsvDataFilePersistence" should be on the server
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BsvDataFilePersistence" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BsvDataFilePersistence gallery card
    Then the "1bdq-obs-pred" table view should open with 22 rows
    And "Reset Camera" button in Biostructure viewer should be visible
    And "Data File" input in Biostructure viewer should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
