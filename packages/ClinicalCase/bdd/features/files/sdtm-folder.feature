Feature: A folder of SDTM data in Browse > Files offers to run Clinical Case
  Claims the package's folder viewer (clinicalCaseFolderLauncher, src/package.ts): a folder holding dm.csv
  shows, on top of its content, that it contains SDTM data, the study id read from dm.csv, and a Run
  ClinicalCase button. The package README has no line on it; the texts are the viewer's own.
  Left out on purpose: the Run ClinicalCase button itself, which opens every domain file as a table view
  and then the application, whose listing writes into App Data (see app/clinical-case-app.feature).

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---ClinicalCase tree node inside browse tree is expanded
    And Files---App-Data---ClinicalCase---SDTM tree node inside browse tree is expanded

  Scenario: The study folder shows the launcher with the study id from its dm.csv
    When user clicks on Files---App-Data---ClinicalCase---SDTM---CDISCPILOT01 tree node inside browse tree
    Then the current view should show the text "Folder contains SDTM data"
    And the current view should show the text "Study ID: CDISCPILOT01"
    And "Run ClinicalCase" button should be visible
    And no errors should have been logged
