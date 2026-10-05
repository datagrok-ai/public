@journey @viewers @realizes:biostructureviewer.viewer.biostructure
Feature: Mol* viewport overlay buttons of the Biostructure viewer
  The row of overlay buttons over the Mol* viewport — Screenshot / State Snapshot, Toggle Controls
  Panel, Toggle Selection Mode, Settings / Controls Info, the package's own Binding site button
  and Toggle Expanded Viewport — each found by its tooltip. Translated from the TestTrack case
  BiostructureViewer/molstar-overlay-extension; Reset Camera is in the smoke feature.

  A button's own on/off state is only a class (`msp-btn-link-toggle-on` / `-off`), which no step
  reads: every "the button class becomes …" line of the md is kept out and requested (see the
  request document). What each button opens is claimed instead, by a control only that panel has:
  Auto-crop for the screenshot panel, Assembly for the structure controls, Mouse Controls for the
  settings panel, "Turn selection mode off" for the selection toolbar. Toggle Expanded Viewport is
  present on the stand (Mol*'s default shows it); its `msp-layout-expanded` class is not readable
  either, so the expanded layout is claimed by the structure controls it brings along, and Escape
  by their going away with the table view still current.

  The Background picks the `pdb` column in Biostructure Id through the settings, because the viewer
  does not take it by itself on the stand (GROK-21119). Adding the
  viewer and building its engine are checked for errors and balloons at the end of the Background.

  Layout Show Controls drives the Mol* layout, but the Toggle Controls Panel button does not write
  back into it: after the button, the property still reads false (checked on the stand, a fact of
  the product, not claimed).

  The Binding site popover is a bare panel attached to the page body, a dialog by its role; its
  header is claimed inside it, not anywhere on the page (the settings have a Binding Site category
  too). The popover stays in the page once opened and is hidden with a class, so after every click
  that opens it its "Show side chains" box is claimed visible before it is read. The box is a bare
  checkbox named by its aria-label, which the library reaches as a text input.

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
    When user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And Biostructure viewer should be visible
    When user clicks on settings icon of Biostructure viewer
    And user selects "pdb" in "Biostructure Id" property in context panel
    Then "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Screenshot / State Snapshot opens its panel and closes it again
    Then "Auto-crop" button in Biostructure viewer should be absent
    When user clicks on "Screenshot / State Snapshot" button in Biostructure viewer
    Then "Auto-crop" button in Biostructure viewer should be visible
    When user clicks on "Screenshot / State Snapshot" button in Biostructure viewer
    Then "Auto-crop" button in Biostructure viewer should be absent
    And no errors should have been logged

  Scenario: Toggle Controls Panel shows and hides the structure controls
    Then "Assembly" button in Biostructure viewer should be absent
    When user clicks on "Toggle Controls Panel" button in Biostructure viewer
    Then "Assembly" button in Biostructure viewer should be visible
    When user clicks on "Toggle Controls Panel" button in Biostructure viewer
    Then "Assembly" button in Biostructure viewer should be absent
    And no errors should have been logged

  Scenario: Layout Show Controls in the settings shows and hides the same controls
    Given "Layout" category in context panel is expanded
    When user checks "Layout Show Controls" property in context panel
    Then "Assembly" button in Biostructure viewer should be visible
    When user unchecks "Layout Show Controls" property in context panel
    Then "Assembly" button in Biostructure viewer should be absent
    And no errors should have been logged

  Scenario: Toggle Selection Mode shows and hides the selection toolbar
    Then "Turn selection mode off" button in Biostructure viewer should be absent
    When user clicks on "Toggle Selection Mode" button in Biostructure viewer
    Then "Turn selection mode off" button in Biostructure viewer should be visible
    When user clicks on "Toggle Selection Mode" button in Biostructure viewer
    Then "Turn selection mode off" button in Biostructure viewer should be absent
    And no errors should have been logged

  Scenario: Settings / Controls Info opens the Mol* settings panel and closes it again
    Then "Mouse Controls" button in Biostructure viewer should be absent
    When user clicks on "Settings / Controls Info" button in Biostructure viewer
    Then "Mouse Controls" button in Biostructure viewer should be visible
    When user clicks on "Settings / Controls Info" button in Biostructure viewer
    Then "Mouse Controls" button in Biostructure viewer should be absent
    And no errors should have been logged

  Scenario: The Binding site button and the Binding Site settings follow each other
    Then "Binding site" button in Biostructure viewer should be enabled
    And "Binding Site" text in dialog should be absent
    When user clicks on "Binding site" button in Biostructure viewer
    Then "Binding Site" text in dialog should be visible
    And "Show side chains" text input should be visible
    And "Show side chains" text input should be unchecked
    And "5.0 Å" text should be visible
    When user checks "Show side chains" text input
    Then "Show side chains" text input should be checked
    And "showBindingSite" property of Biostructure viewer should be "true"
    Given "Binding Site" category in context panel is expanded
    Then "Show Binding Site" property in context panel should be checked
    When user enters "8" in "Binding Site Radius" property in context panel
    Then "bindingSiteRadius" property of Biostructure viewer should be "8"
    When user clicks on "Binding site" button in Biostructure viewer
    Then "Show side chains" text input should be visible
    And "8.0 Å" text should be visible
    When user unchecks "Show Binding Site" property in context panel
    Then "showBindingSite" property of Biostructure viewer should be "false"
    When user clicks on "Binding site" button in Biostructure viewer
    Then "Show side chains" text input should be visible
    And "Show side chains" text input should be unchecked
    When user presses Escape
    Then "Show side chains" text input should be hidden
    When user enters "5" in "Binding Site Radius" property in context panel
    Then "bindingSiteRadius" property of Biostructure viewer should be "5"
    And no errors should have been logged

  Scenario: Toggle Expanded Viewport expands the Mol* layout and Escape brings it back
    Then "Assembly" button in Biostructure viewer should be absent
    When user clicks on "Toggle Expanded Viewport" button in Biostructure viewer
    Then "Assembly" button in Biostructure viewer should be visible
    When user presses Escape
    Then "Assembly" button in Biostructure viewer should be absent
    And the "pdb_data" view should be current
    And Biostructure viewer should be visible
    And no errors should have been logged
