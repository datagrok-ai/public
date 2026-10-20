@journey @realizes:views.projects
Feature: A project saved with and without Data sync
  demog.csv opened from Browse > Files > Demo is saved through the ribbon's Save dialog with Data
  sync on, then saved as a copy with Data sync off, and the copy is saved again with Data sync on.
  What the server keeps is claimed for every save (the table's data-sync flag and creation script).
  Translated from the TestTrack case Projects/complex-save-copy.

  How a reopened project brings the table back — rebuilt by its creation script, or loaded from the
  uploaded snapshot — is ApiTests `Dapi: projects: data sync` (src/dapi/projects.ts).

  Both projects are named with the run's time and removed, with their tables and views, when the
  feature ends.

  Not translated, and why: opening the projects from the Dashboards tiles (the tile gallery is
  Browse's subject) and the Share dialog the first save
  shows, which is cancelled (a copy of a saved project shows none).

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDD-Sync-Orig-{time}" is on the server
    And no project named "BDD-Sync-Copy-{time}" is on the server

  Scenario: A file opened from Browse is saved with Data sync on
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 5850 rows
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDD-Sync-Orig-{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDD-Sync-Orig-{time}" should have been shown
    And 1 project named "BDD-Sync-Orig-{time}" should be on the server
    And the "demog" table of the "BDD-Sync-Orig-{time}" project should be saved with data sync
    And "Share BDD-Sync-Orig-{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDD-Sync-Orig-{time}" dialog
    Then the "Share BDD-Sync-Orig-{time}" dialog should close
    And no errors should have been logged

  Scenario: A copy saved with Data sync off leaves the original synced
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of BDD-Sync-Orig-{time}"
    And choice input in "demog" project table in "Save project" dialog should have value "Clone"
    When user enters "BDD-Sync-Copy-{time}" into Name text input in "Save project" dialog
    And user unchecks Data sync switch in "demog" project table in "Save project" dialog
    Then Data sync switch in "demog" project table in "Save project" dialog should be unchecked
    And "Creation script" button in "demog" project table in "Save project" dialog should be hidden
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDD-Sync-Copy-{time}" should be on the server
    And the "demog" table of the "BDD-Sync-Copy-{time}" project should be saved as a snapshot
    And the "demog" table of the "BDD-Sync-Orig-{time}" project should be saved with data sync
    And no errors should have been logged

  Scenario: The copy's Save dialog keeps Data sync off
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be unchecked
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged

  Scenario: Turning Data sync on for the copy saves its creation script
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user checks Data sync switch in "demog" project table in "Save project" dialog
    And user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then "demog" project table in "Save project" dialog should contain text "System:DemoFiles/demog.csv"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "demog" table of the "BDD-Sync-Copy-{time}" project should be saved with data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown
