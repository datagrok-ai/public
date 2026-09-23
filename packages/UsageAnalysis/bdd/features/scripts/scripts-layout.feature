@journey @serial @realizes:views.scripts
Feature: The layout of a script's result
  A script that outputs a dataframe gets a Layout tab: run it there, arrange the viewers the result
  opens with, and Save stores them with the script. Translated from files/TestTrack/Scripts/layout.md
  and playwright-public/scripts/scripts-layout.test.ts, which reapplied the layout through the JS
  API and skipped its own save when the dialog did not close — so it proved the API, not the view.

  The script is this feature's own ({time} in its name) and is deleted at the end, together with the
  layout the save writes.

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

  Scenario: Save stores the layout with the script
    When user saves the script
    Then no error or warning balloon should have been shown
    And the script "BddScriptLayout{time}" on the server should have a layout
    And no errors should have been logged
