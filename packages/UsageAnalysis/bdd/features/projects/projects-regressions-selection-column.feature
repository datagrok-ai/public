@realizes:views.projects
Feature: Projects regressions: a column saved from the row selection
  core#602: demog.csv from Browse > Files > Demo gets ten rows selected, and Selected > Save as Column
  on the status bar adds the bool column A that marks them. The project is saved through the ribbon's
  Save dialog and reopened from its Dashboards card. Saved as a snapshot, the table comes back with A.
  Saved with Data sync, the table is rebuilt by its creation script, which does not record the
  column, and it comes back without A; the Save dialog says nothing about it. That scenario is a known
  failure until the column survives the reopen or the dialog warns before the save.

  The selection is set through the API: the rows are a precondition, the gesture under test is Save
  as Column. The project is named with the run's time and removed with its table and view when each
  scenario starts and when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user selects the first 10 rows
    And user clicks on "Selected: 10" text in status bar
    And user picks "Save as Column" from the open menu
    Then "A" column should have type "bool"

  Scenario: Saved as a snapshot, the table comes back with the selection column
    Given no project named "BDDRegSelColSnap{time}" is on the server
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegSelColSnap{time}" into Name text input in "Save project" dialog
    And user unchecks Data sync switch in "demog" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegSelColSnap{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegSelColSnap{time}" dialog
    Then the "Share BDDRegSelColSnap{time}" dialog should close
    And the "demog" table of the "BDDRegSelColSnap{time}" project should be saved as a snapshot
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegSelColSnap{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegSelColSnap{time} gallery card and the project opens
    Then the "demog" view should be current
    And the table should have been loaded as a snapshot
    And "A" column should have type "bool"

  # core#602: the creation script does not record Save as Column, and the dialog does not warn
  @known-failure:70 @realizes:core#602
  Scenario: Saved with Data sync, the table comes back with the selection column
    Given no project named "BDDRegSelColSync{time}" is on the server
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDRegSelColSync{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegSelColSync{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegSelColSync{time}" dialog
    Then the "Share BDDRegSelColSync{time}" dialog should close
    And the "demog" table of the "BDDRegSelColSync{time}" project should be saved with data sync
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegSelColSync{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegSelColSync{time} gallery card and the project opens
    Then the "demog" view should be current
    And the table should have been reloaded by data sync
    And "A" column should have type "bool"
