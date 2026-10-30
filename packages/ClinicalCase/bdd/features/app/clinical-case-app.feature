@journey
Feature: Clinical Case opens from Browse > Apps with its demo studies and the Import study view
  Claims the start of the application as the package README describes it ("Initial application page
  contains the list of pre-loaded studies"; "Click Import study tree item"): the start page with the two
  demo SDTM studies of the package's files (CDISCPILOT01, nida-ctn-0001) and their subject counts, the
  studies and Import study under the Clinical Case node of the tree, and the Import study view with its
  Import button held disabled, and the error it shows for files without a study definition.
  Listing the studies writes a studyConfig.json into each study folder of System:AppData/ClinicalCase; the
  first scenario, right after the listing, puts the same text there through the library's file step, which
  deletes both files at feature end (a run killed before that leaves the same text the app writes).
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
    # the listing above wrote these two files; the step deletes them at feature end
    Given a file "System:AppData/ClinicalCase/SDTM/CDISCPILOT01/studyConfig.json" with text '{"standard":"SDTM","name":"CDISCPILOT01","totalSubjects":306}' is on the server
    And a file "System:AppData/ClinicalCase/SDTM/nida-ctn-0001/studyConfig.json" with text '{"standard":"SDTM","name":"nida-ctn-0001","totalSubjects":138,"other":{"TEST PRODUCT IS BEING ADDED ON TO EXISTING TREATMENT":"Y","AGE GROUP":"CHILDREN AND ADULTS","PLANNED MAXIMUM AGE OF SUBJECTS":"NO UPPER LIMIT","PLANNED MINIMUM AGE OF SUBJECTS":"15 YEARS","COMPARATIVE TREATMENT NAME":"CLONIDINE","TYPE OF CONTROL":"ACTIVE","TRIAL DESIGN":"SINGLE GROUP PARALLEL","DIAGNOSIS GROUP":"ADDICTED SUBJECTS","TEST PRODUCT TOTAL DAILY DOSE":"ESCALATING FROM 4 MG BUP/1 MG NX TO 16 MG BUP TO 4 MG NX","TEST PRODUCT DOSE UNITS":"MG","DRUG DICTIONARY":"NONE","TRIAL INDICATIONS":"OPIATE DEPENDENCE","TRIAL INDICATION TYPE":"TREATMENT","LENGTH OF TRIAL":"P6M","TRIAL BLINDING SCHEMA":"OPEN LABEL","TRIAL PRIMARY OBJECTIVE":"COMPARE THE RELATIVE CLINICAL UTILITY OF BUP/NX TO CLONIDINE IN A SHORT-TERM (13 DAY) DETOXIFICATION INITIATED IN AN INPATIENT SETTING","TRIAL SECONDARY OBJECTIVE":"DETERMINE THE EFFECT OF TREATMENT ON WITHDRAWAL SYMPTOMS, TREATMENT RETENTION, SAFETY, PATIENT SATISFACTION, AND EXTENT OF ANCILLARY MEDICATIONS DISPENSED","TRIAL PHASE":"4","PLANNED NUMBER OF SUBJECTS":"360","TRIAL IS RANDOMIZED":"Y","RANDOMIZATION RATIO":"2:1 TEST TO COMPARATOR","TEST PRODUCT ROUTE OF ADMINISTRATION":"SUBLINGUAL","SEX OF PARTICIPANTS":"BOTH","SPONSORING ORGANIZATION":"NIDA","TRIAL TITLE":"BUPRENORPHINE/NALOXONE VERSUS CLONIDINE FOR INPATIENT OPIATE DETOXIFICATION","REPORTED NAME OF TEST PRODUCT":"BUPRENORPHINE/CLONIDINE","TYPE OF TRIAL":"EFFICACY"}}' is on the server
    Then no errors should have been logged

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
