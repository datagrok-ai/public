@journey @serial @realizes:views.scripts
Feature: The layout of a script's result
  A script that outputs a dataframe gets a Layout tab: run it there, arrange the viewers, colours and
  columns the result opens with, and Save stores them with the script; a run of the saved script opens
  its result that way. Translated from files/TestTrack/Scripts/layout.md
  and playwright-public/scripts/scripts-layout.test.ts, which reapplied the layout through the JS
  API and skipped its own save when the dialog did not close — so it proved the API, not the view.

  The result of a run, with one more viewer, is then saved as a project through the ribbon's Save
  dialog and opened again with its viewers, hidden column and colouring, which the Source pane's
  REFRESH (the md's "Toolbox > File > Refresh") keeps (md steps 9-14). REFRESH is shown to run the
  script again by a deleted row coming back.

  The script and the project are this feature's own ({time} in their names) and are deleted at the
  end, together with the layout the save writes.

  Not translated, and why: docking one viewer over another (md step 4) has no named drop target —
  the dock manager's drop zones carry no name, so there is nothing for a gesture to aim at.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And a script "BddScriptLayout{time}" is on the server:
      """
      //language: javascript
      //output: dataframe df
      df = await grok.data.getDemoTable('cars.csv');
      """
    And the layouts saved for the script are deleted at the end
    And no project named "BDD-Script-layout-proj-{time}" is on the server
    And user opens the Scripts view

  Scenario: The Layout tab asks for a run first
    When user types "BddScriptLayout{time}" into gallery search
    And user double-clicks on "BddScriptLayout{time}" link in gallery
    Then the "BddScriptLayout{time}" view should be current
    When user clicks on Layout tab
    Then "Run script to get data and edit layout." text should be visible
    And grid viewer should be absent
    And no errors should have been logged

  Scenario: Running the script fills the layout with its result
    When user clicks on "Run script (F5)" icon
    Then grid viewer should be visible
    And the "rows" reading of grid should be 30
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A viewer added from the toolbox joins the layout
    When user clicks on "bar chart" icon in toolbox
    Then bar chart viewer should be visible
    And bar chart viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A colouring and a hidden column join the layout
    Then the "color of cell 1 of length" and "color of cell 2 of length" readings of grid should be the same
    When user picks "Color Coding > Linear" from the context menu of the "header length" area of grid
    Then the "color of cell 1 of length" and "color of cell 2 of length" readings of grid should differ
    When user picks "Hide" from the context menu of the "header model" area of grid
    Then grid should not have a "header model" area
    And no errors should have been logged

  Scenario: Save stores the layout with the script
    When user saves the script
    Then no error or warning balloon should have been shown
    And the script "BddScriptLayout{time}" on the server should have a layout
    And no errors should have been logged

  Scenario: Run... from the Scripts view opens the result with the saved layout
    When user closes all views
    Given user opens the Scripts view
    When user types "BddScriptLayout{time}" into gallery search
    And user picks "Run..." from the context menu of "BddScriptLayout{time}" link in gallery
    Then the current view should be a TableView view
    And the table should have 30 rows
    And the open tableview should have 1 bar chart viewer
    And grid should not have a "header model" area
    And the "color of cell 1 of length" and "color of cell 2 of length" readings of grid should differ
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The result with one more viewer is saved as a project and opens with its layout
    Given the toolbox pane is shown
    When user clicks on histogram icon in toolbox
    Then the open tableview should have 1 histogram viewer
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user enters "BDD-Script-layout-proj-{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDD-Script-layout-proj-{time}" should be on the server
    And "Share BDD-Script-layout-proj-{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDD-Script-layout-proj-{time}" dialog
    Then the "Share BDD-Script-layout-proj-{time}" dialog should close
    When user closes all views
    And user opens the "BDD-Script-layout-proj-{time}" project and waits for its table
    Then the table should have 30 rows
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 histogram viewer
    And grid should not have a "header model" area
    And the "color of cell 1 of length" and "color of cell 2 of length" readings of grid should differ
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: REFRESH of the reopened result runs the script again and keeps its layout
    # a deleted row is what shows the script ran again: REFRESH brings it back
    When user clicks on the "row header 1" area of grid holding Shift
    And user presses Shift+Delete
    Then the table should have 29 rows
    Given the toolbox pane is shown
    When user clicks on REFRESH button in toolbox
    Then the table should have 30 rows
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 histogram viewer
    And grid should not have a "header model" area
    And the "color of cell 1 of length" and "color of cell 2 of length" readings of grid should differ
    And no errors should have been logged
    And no error or warning balloon should have been shown
