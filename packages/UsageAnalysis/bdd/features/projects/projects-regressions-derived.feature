@serial @realizes:views.projects
Feature: Projects regressions: a data-sync project with a published pivot saved as a zip
  Regression guard for GROK-19580 (Save as Zip of a data-sync project with a published pivot): demog
  is opened by a double-click in Browse > Files > Demo (it stands for SPGI in the report), a pivot
  table from the toolbox is published into the workspace, the project is saved with Data sync
  through the ribbon's Save dialog and saved as a zip from its Dashboards card. The scenario fails
  on no zip downloaded, a zip without the pivot table, or an error balloon.

  Parked in the request document until the library has the phrases: GROK-19664 (a pivot
  published from a second table of the same name goes into the other open project), which reads the
  Dashboards panel of the left sidebar and the tables open in the workspace, and GROK-19135 (linking
  works only in the first of two open projects with the same linked tables), which reads a table of
  one open project among two that hold tables of the same names. GROK-20013 (a copy saved after an
  in-place join) is parked with projects-regressions-join.

  The console errors are not claimed: the Save dialog's preview logs "Unable to find element in
  cloned iframe" for these views on every run (GROK-18606, won't fix), and the check that lets only
  that message through is requested in the request document. The project is named with the run's time and
  removed with its tables and views when the scenario starts and when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-19580
  Scenario: A data-sync project with a published pivot is saved as a zip from the gallery
    Given no project named "BDDRegZip{time}" is on the server
    And user clears the saved pivot table parameters
    Given the browse panel is open
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog aggregation" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    And Data sync switch in "demog aggregation" project table in "Save project" dialog should be checked
    When user enters "BDDRegZip{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegZip{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegZip{time}" dialog
    Then the "Share BDDRegZip{time}" dialog should close
    And the "demog aggregation" table of the "BDDRegZip{time}" project should be saved with data sync
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegZip{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches downloads
    When user picks "Save as Zip" from the context menu of BDDRegZip{time} gallery card
    Then a file "BDDRegZip{time}.zip" should have been downloaded
    And the downloaded file "BDDRegZip{time}.zip" should contain the text "DemogAggregation = Aggregate("
    And no error or warning balloon should have been shown
