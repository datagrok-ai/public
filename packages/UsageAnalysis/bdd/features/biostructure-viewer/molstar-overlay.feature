@journey @viewers @realizes:biostructureviewer.viewer.biostructure
Feature: Mol* viewport overlay buttons of the Biostructure viewer
  The row of overlay buttons over the Mol* viewport — Screenshot / State Snapshot, Toggle Controls
  Panel, Toggle Selection Mode, Settings / Controls Info, the package's own Binding site button
  and Toggle Expanded Viewport — each found by its tooltip. Translated from the TestTrack case
  BiostructureViewer/molstar-overlay-extension; Reset Camera is in the smoke feature.

  A button's own on/off state is Mol*'s class `msp-btn-link-toggle-on`, which the library reads as
  selected, and what each button opens is claimed as well, by a control only that panel has:
  Auto-crop for the screenshot panel, Assembly for the structure controls, Mouse Controls for the
  settings panel, "Turn selection mode off" for the selection toolbar. The expanded layout is the
  viewer's `layout expanded` reading, and Escape brings it back with the table view still current.

  The Background picks the `pdb` column in Biostructure Id through the settings, because the viewer
  does not take it by itself on the stand (GROK-21119). Adding the
  viewer and building its engine are checked for errors and balloons at the end of the Background.

  Layout Show Controls drives the Mol* layout, but the Toggle Controls Panel button does not write
  back into it: after the button, the property still reads false (checked on the stand, a fact of
  the product, not claimed).

  The Binding site popover is a bare panel attached to the page body, a dialog by its role; its
  header is claimed inside it, not anywhere on the page (the settings have a Binding Site category
  too). The popover stays in the page once opened and is hidden with a class, so after every click
  that opens it its "Show side chains" box (a bare checkbox named by its aria-label) is claimed
  visible before it is read. A new radius is claimed on the atoms of the binding site Mol* built
  (335 at 5 Å, 559 at 8 Å), not on the property the step just wrote.

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
    When user types "Biostructure" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And Biostructure viewer should be visible
    When user clicks on settings icon of Biostructure viewer
    And user selects "pdb" in "Biostructure Id" property in context panel
    Then "Reset Camera" button in Biostructure viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Screenshot / State Snapshot opens its panel and closes it again
    Then "Auto-crop" button in Biostructure viewer should be absent
    And "Screenshot / State Snapshot" overlay button in Biostructure viewer should not be selected
    When user clicks on "Screenshot / State Snapshot" overlay button in Biostructure viewer
    Then "Auto-crop" button in Biostructure viewer should be visible
    And "Screenshot / State Snapshot" overlay button in Biostructure viewer should be selected
    When user clicks on "Screenshot / State Snapshot" overlay button in Biostructure viewer
    Then "Auto-crop" button in Biostructure viewer should be absent
    And "Screenshot / State Snapshot" overlay button in Biostructure viewer should not be selected
    And no errors should have been logged

  Scenario: Toggle Controls Panel shows and hides the structure controls
    Then "Assembly" button in Biostructure viewer should be absent
    And "Toggle Controls Panel" overlay button in Biostructure viewer should not be selected
    When user clicks on "Toggle Controls Panel" overlay button in Biostructure viewer
    Then "Assembly" button in Biostructure viewer should be visible
    And "Toggle Controls Panel" overlay button in Biostructure viewer should be selected
    And the "controls shown" reading of Biostructure viewer should be "true"
    When user clicks on "Toggle Controls Panel" overlay button in Biostructure viewer
    Then "Assembly" button in Biostructure viewer should be absent
    And "Toggle Controls Panel" overlay button in Biostructure viewer should not be selected
    And no errors should have been logged

  Scenario: Layout Show Controls in the settings shows and hides the same controls
    Given "Layout" category in context panel is expanded
    When user checks "Layout Show Controls" property in context panel
    Then "Assembly" button in Biostructure viewer should be visible
    And "Toggle Controls Panel" overlay button in Biostructure viewer should be selected
    When user unchecks "Layout Show Controls" property in context panel
    Then "Assembly" button in Biostructure viewer should be absent
    And "Toggle Controls Panel" overlay button in Biostructure viewer should not be selected
    And no errors should have been logged

  Scenario: Toggle Selection Mode shows and hides the selection toolbar
    Then "Turn selection mode off" button in Biostructure viewer should be absent
    When user clicks on "Toggle Selection Mode" overlay button in Biostructure viewer
    Then "Turn selection mode off" button in Biostructure viewer should be visible
    And "Toggle Selection Mode" overlay button in Biostructure viewer should be selected
    When user clicks on "Toggle Selection Mode" overlay button in Biostructure viewer
    Then "Turn selection mode off" button in Biostructure viewer should be absent
    And "Toggle Selection Mode" overlay button in Biostructure viewer should not be selected
    And no errors should have been logged

  Scenario: Settings / Controls Info opens the Mol* settings panel and closes it again
    Then "Mouse Controls" button in Biostructure viewer should be absent
    When user clicks on "Settings / Controls Info" overlay button in Biostructure viewer
    Then "Mouse Controls" button in Biostructure viewer should be visible
    And "Settings / Controls Info" overlay button in Biostructure viewer should be selected
    When user clicks on "Settings / Controls Info" overlay button in Biostructure viewer
    Then "Mouse Controls" button in Biostructure viewer should be absent
    And "Settings / Controls Info" overlay button in Biostructure viewer should not be selected
    And no errors should have been logged

  Scenario: The Binding site button and the Binding Site settings follow each other
    Then "Binding site" overlay button in Biostructure viewer should be enabled
    And "Binding site" overlay button in Biostructure viewer should not be selected
    And "Binding Site" text in dialog should be absent
    When user clicks on "Binding site" overlay button in Biostructure viewer
    Then "Binding Site" text in dialog should be visible
    And "Show side chains" checkbox should be visible
    And "Show side chains" checkbox should be unchecked
    And "5.0 Å" text should be visible
    When user checks "Show side chains" checkbox
    Then "Show side chains" checkbox should be checked
    And "showBindingSite" property of Biostructure viewer should be "true"
    And "Binding site" overlay button in Biostructure viewer should be selected
    And the "binding site shown" reading of Biostructure viewer should be "true"
    And the "binding site atoms" reading of Biostructure viewer should be 335
    Given "Binding Site" category in context panel is expanded
    Then "Show Binding Site" property in context panel should be checked
    When user enters "8" in "Binding Site Radius" property in context panel
    Then "bindingSiteRadius" property of Biostructure viewer should be "8"
    And the "binding site atoms" reading of Biostructure viewer should be 559
    When user clicks on "Binding site" overlay button in Biostructure viewer
    Then "Show side chains" checkbox should be visible
    And "8.0 Å" text should be visible
    When user unchecks "Show Binding Site" property in context panel
    Then "showBindingSite" property of Biostructure viewer should be "false"
    And "Binding site" overlay button in Biostructure viewer should not be selected
    And the "binding site shown" reading of Biostructure viewer should be "false"
    When user clicks on "Binding site" overlay button in Biostructure viewer
    Then "Show side chains" checkbox should be visible
    And "Show side chains" checkbox should be unchecked
    When user presses Escape
    Then "Show side chains" checkbox should be hidden
    When user enters "5" in "Binding Site Radius" property in context panel
    Then "bindingSiteRadius" property of Biostructure viewer should be "5"
    And no errors should have been logged

  Scenario: Toggle Expanded Viewport expands the Mol* layout and Escape brings it back
    Then "Assembly" button in Biostructure viewer should be absent
    And the "layout expanded" reading of Biostructure viewer should be "false"
    When user clicks on "Toggle Expanded Viewport" overlay button in Biostructure viewer
    Then the "layout expanded" reading of Biostructure viewer should be "true"
    And "Toggle Expanded Viewport" overlay button in Biostructure viewer should be selected
    And "Assembly" button in Biostructure viewer should be visible
    When user presses Escape
    Then the "layout expanded" reading of Biostructure viewer should be "false"
    And "Assembly" button in Biostructure viewer should be absent
    And the "pdb_data" view should be current
    And Biostructure viewer should be visible
    And no errors should have been logged
