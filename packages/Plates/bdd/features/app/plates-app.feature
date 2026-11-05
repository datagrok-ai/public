@journey
Feature: Apps > Plates: the plate list, the searches, Templates and Create Plate open with their content
  Apps > Plates (`platesApp`, `initPlatesAppTree` in `src/plates/plates-app.ts`) opens a table of the
  registered plates, one column per plate property, and fills its tree node with Create, Search
  plates, Search wells, Search analyses and Templates. There is no help page for Plates; what is
  claimed is the package's own behaviour. The plates are the demo set the package seeds into its own
  database on first use (`plates-demo.ts`): 20 plates of the "Cell counting" template and 20 of
  "Dose-response", 40 in all — the counts below are that set's, measured on the local stand, and hold
  only while the demo set is the only data in the package database (a plate created by hand or by
  another test changes them).

  Left out on purpose:
  - Search Plates and Search Analyses opening with their inputs: `UsageAnalysis` browse-apps-matrix
    claims them; here only the search itself (a number typed into a property narrows the plates);
  - SAVE in a template editor, CREATE with an imported plate, a template's Clone / Delete: each writes
    plates or templates into the package database, which no library step puts back (MISSING.md);
  - the values a search finds by a choice property: the demo plates get random values, so no exact
    count exists; the matchers are `src/plates/matchers.ts`'s job (package tests);
  - the Assay Plates demo: its Excel plates (`System:DemoFiles/hts/xlsx_plates`) are not on the local
    stand (MISSING.md).
  The scenarios save nothing themselves, but the first open of the app seeds the two demo templates and
  their 40 plates into the package database when they are not there yet (`seedDemoPlateData`). Each
  scenario opens its view from the tree again: the app's own view hides the Browse panel, and a click on
  the node already selected opens nothing.

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Plates tree node inside browse tree is expanded

  Scenario: Apps > Plates lists the registered plates with their plate properties as columns
    When user clicks on Apps---Plates tree node inside browse tree
    Then the "Plates" view should be current
    And the table should have 40 rows
    And "barcode" column should have semantic type "Plate Barcode"
    And the table should have a column "Imaging device"
    And the table should have a column "Project"
    And no errors should have been logged

  Scenario: Search Wells offers the well properties of the template picked in it
    Given the browse panel is open
    And Apps---Plates tree node inside browse tree is expanded
    When user clicks on Apps---Plates---Search-wells tree node inside browse tree
    Then the "Search Wells" view should be current
    And "Template" input should have value "Cell counting"
    And "Well cell count" input should be visible
    When user selects "Dose-response" in "Template" input
    Then "Concentration" input should be visible
    And "Activity" input should be visible
    And "Well cell count" input should be hidden
    When user selects "Cell counting" in "Template" input
    Then "Well cell count" input should be visible
    And "Concentration" input should be hidden
    And no errors should have been logged

  Scenario: Search Plates narrows the plates found by a number typed into a plate property
    Given the browse panel is open
    And Apps---Plates tree node inside browse tree is expanded
    When user clicks on Apps---Plates---Search-plates tree node inside browse tree
    Then the "Search Plates" view should be current
    And the "rows" reading of grid should be 40
    When user enters ">=0" into "Plate cell count" input
    Then the "rows" reading of grid should be 20
    When user clears "Plate cell count" input
    Then the "rows" reading of grid should be 40
    And no errors should have been logged

  Scenario: Templates lists the templates; EDIT opens one, CREATE NEW a new one to save
    Given the browse panel is open
    And Apps---Plates tree node inside browse tree is expanded
    When user clicks on Apps---Plates---Templates tree node inside browse tree
    Then the "Templates" view should be current
    And the current view should show the text "Cell counting"
    And the current view should show the text "Dose-response"
    When user clicks on first "EDIT" button
    Then the "Templates / Cell counting" view should be current
    And "Name" input should have value "Cell counting"
    And "SAVE" button should be hidden
    When user clicks on Apps---Plates---Search-wells tree node inside browse tree
    And user clicks on Apps---Plates---Templates tree node inside browse tree
    Then the "Templates" view should be current
    When user clicks on "CREATE NEW" button
    Then the "Templates / New Template" view should be current
    And "Name" input should have value "New Template"
    And "SAVE" button should be visible
    And no errors should have been logged

  Scenario: Create Plate shows an empty plate with an analysis tab each, and CREATE with no plate imported only warns
    Given the browse panel is open
    And Apps---Plates tree node inside browse tree is expanded
    When user clicks on Apps---Plates---Create tree node inside browse tree
    Then the "Create Plate" view should be current
    And "Plate View" tab should be visible
    And "qPCR Analysis" tab should be visible
    When user clicks on "Dose Response" tab
    Then the current view should show the text "To see dose response curves, please map the following required properties: Activity, Concentration, SampleID."
    When user clicks on "CREATE" button
    Then a warning balloon containing "No active plate is selected." should have been shown
    And no errors should have been logged
