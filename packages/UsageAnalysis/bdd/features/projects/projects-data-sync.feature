@journey @serial @realizes:views.projects
Feature: A project saved with and without Data sync
  demog.csv opened from Browse > Files > Demo is saved through the ribbon's Save dialog with Data
  sync on, then saved as a copy with Data sync off; the Dashboards panel of the left sidebar names
  the copy as the open project. The copy reopens from its Dashboards card with no creation script,
  is saved again with Data sync on and then shows the file's OpenFile call as its creation script,
  while the original keeps its own. What the server keeps is claimed for every save (the table's
  data-sync flag), and every reopen by how the table came back: rebuilt by its creation script, or
  loaded from the uploaded snapshot. Translated from the TestTrack case Projects/complex-save-copy.

  Nothing of the md is parked. Names are letters and digits only (the Dashboards search misses "-"
  and "_") and carry the run's time; both projects are removed, with their tables and views, when
  the feature starts and ends. It is @serial: a save uploads the table and the reopen reads it back.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDSyncOrig{time}" is on the server
    And no project named "BDDSyncCopy{time}" is on the server

  Scenario: A file opened from Browse is saved with Data sync on
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 5850 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDSyncOrig{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDSyncOrig{time}" uploaded' should have been shown
    And 1 project named "BDDSyncOrig{time}" should be on the server
    And the "demog" table of the "BDDSyncOrig{time}" project should be saved with data sync
    And "Share BDDSyncOrig{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDSyncOrig{time}" dialog
    Then the "Share BDDSyncOrig{time}" dialog should close
    And no errors should have been logged

  Scenario: A copy saved with Data sync off leaves the original synced, and is the open project
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of BDDSyncOrig{time}"
    And choice input in "demog" project table in "Save project" dialog should have value "Clone"
    When user enters "BDDSyncCopy{time}" into Name text input in "Save project" dialog
    And user unchecks Data sync switch in "demog" project table in "Save project" dialog
    Then "Creation script" button in "demog" project table in "Save project" dialog should be hidden
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDSyncCopy{time}" uploaded' should have been shown
    And 1 project named "BDDSyncCopy{time}" should be on the server
    And the "demog" table of the "BDDSyncCopy{time}" project should be saved as a snapshot
    And the "demog" table of the "BDDSyncOrig{time}" project should be saved with data sync
    When user clicks on Dashboards tab
    Then BDDSyncCopy{time} tree node inside browse tree should be visible
    And BDDSyncOrig{time} tree node inside browse tree should be absent
    # the sidebar tab toggles the Dashboards panel; left open, it adds a second Save button to every table view
    When user clicks on Dashboards tab
    Then "New Dashboard" tree node inside browse tree should be hidden
    And no errors should have been logged

  Scenario: The copy opens from its card as a snapshot, with no creation script
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSync" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then BDDSyncOrig{time} gallery card should be visible
    And BDDSyncCopy{time} gallery card should be visible
    When user double-clicks on BDDSyncCopy{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been loaded as a snapshot
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be unchecked
    And "Creation script" button in "demog" project table in "Save project" dialog should be hidden
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged

  Scenario: Data sync turned on for the copy, the copy rereads the file
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user checks Data sync switch in "demog" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "demog" table of the "BDDSyncCopy{time}" project should be saved with data sync
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSyncCopy{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDSyncCopy{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then "demog" project table in "Save project" dialog should contain text 'OpenFile("System:DemoFiles/demog.csv")'
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged

  Scenario: The original still opens by rereading the file, with its creation script
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSyncOrig{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDSyncOrig{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown
