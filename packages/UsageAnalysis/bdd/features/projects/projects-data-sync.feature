@journey @serial @realizes:views.projects
Feature: A project saved with and without Data sync
  demog.csv opened from Browse > Files > Demo is saved through the ribbon's Save dialog with Data
  sync on, then saved as a copy with Data sync off, and the copy is saved again with Data sync on.
  What the server keeps is claimed for every save (the table's data-sync flag and creation script),
  and every reopen is claimed by how the table came back: rebuilt by its creation script, or loaded
  from the uploaded snapshot. Translated from the TestTrack case Projects/complex-save-copy.

  Both projects are named with the run's time and removed, with their tables and views, when the
  feature ends. It is @serial: a save uploads the table and the reopen reads it back.

  Not translated, and why: opening the projects from the Dashboards tiles (the tile gallery is
  Browse's subject; the projects open through the API here) and the Share dialog the first save
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
    Then "Creation script" button in "demog" project table in "Save project" dialog should be hidden
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDD-Sync-Copy-{time}" should be on the server
    And the "demog" table of the "BDD-Sync-Copy-{time}" project should be saved as a snapshot
    And the "demog" table of the "BDD-Sync-Orig-{time}" project should be saved with data sync
    And no errors should have been logged

  Scenario: The copy opens from its snapshot
    When user closes all views
    And user opens the "BDD-Sync-Copy-{time}" project and waits for its table
    Then the table should have 5850 rows
    And the table should have been loaded as a snapshot
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be unchecked
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged

  Scenario: Turning Data sync on for the copy makes it reread the file
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user checks Data sync switch in "demog" project table in "Save project" dialog
    And user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then "demog" project table in "Save project" dialog should contain text "System:DemoFiles/demog.csv"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "demog" table of the "BDD-Sync-Copy-{time}" project should be saved with data sync
    When user closes all views
    And user opens the "BDD-Sync-Copy-{time}" project and waits for its table
    Then the table should have 5850 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged

  Scenario: The original still opens by rereading the file
    When user closes all views
    And user opens the "BDD-Sync-Orig-{time}" project and waits for its table
    Then the table should have 5850 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown
