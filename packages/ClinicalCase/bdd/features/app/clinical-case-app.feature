@journey
Feature: Clinical Case opens from Browse > Apps with its demo studies and the Import study view
  Claims the start of the application as the package README describes it ("Initial application page
  contains the list of pre-loaded studies"; "Click Import study tree item"): the start page with the two
  demo SDTM studies of the package's files (CDISCPILOT01, nida-ctn-0001) and their subject counts, the
  studies and Import study under the Clinical Case node of the tree, and the Import study view with its
  Import button held disabled, and the error it shows for files without a study definition.
  Listing the studies writes a studyConfig.json, the app's own cache, into each demo study folder of
  System:AppData/ClinicalCase; the feature leaves it there, as a user's first visit does.
  Left out on purpose: opening a study (its Summary and every view under it) and Import itself. Opening a
  study runs the CDISC rules engine in the package's container (ClinicalCase:RunCoreValidate) and writes
  <study>.d42 and validation_results.json into the study folder, which no library step puts back; Import
  writes the chosen files there. See MISSING.md.

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded

  Scenario: The start page lists the demo studies with their subject counts
    When user clicks on Apps---Clinical-Case tree node inside browse tree
    Then the "Clinical Case" view should be current
    And the current view should show the text "Visualize and explore your SDTM data"
    And the current view should show the text "Study Standard Protocol Total subjects Start date End date CDISCPILOT01 SDTM 306 nida-ctn-0001 SDTM 138"
    And "Import study..." button should be visible
    And no errors should have been logged

  Scenario: The Clinical Case node lists the studies and Import study
    When user expands Apps---Clinical-Case tree node inside browse tree
    Then Apps---Clinical-Case---CDISCPILOT01 tree node inside browse tree should be visible
    And Apps---Clinical-Case---nida-ctn-0001 tree node inside browse tree should be visible
    And Apps---Clinical-Case---Import-study tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Import study... opens the import view, which refuses files without a study definition
    When user clicks on Apps---Clinical-Case tree node inside browse tree
    Then the "Clinical Case" view should be current
    When user clicks on "Import study..." button
    Then the "Import Study - Clinical Case" view should be current
    And the current view should show the text "Import study files"
    And "Import" button should be disabled
    When user uploads "../files/SDTM/CDISCPILOT01/dm.csv" through "Open files" icon
    Then the current view should show the text "Error: Invalid or missing define.xml/study.json"
    And the current view should show the text "dm.csv"
    And an error balloon containing "Invalid or missing define.xml/study.json" should have been shown
    And "Import" button should be disabled
    When user closes the current view
    Then no errors should have been logged
