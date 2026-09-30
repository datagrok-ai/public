@journey @serial @realizes:views.projects @realizes:sharing.share-dialog
Feature: A project built from a file, shared, renamed and deleted
  demog.csv opened from Browse > Files > Demo is saved as a project with Data sync through the
  ribbon's Save dialog; the project is shared with the second account (View and use, notifications
  off) from its Dashboards card, renamed from the card and reopened under its new name, and deleted
  through Delete Project. Translated from the TestTrack case Projects/projects-lifecycle-files.

  Parked (see the request document): the recipient's opens (with View and use it cannot save the
  original; with Full access it can), which need signing in as the second account, and the grant
  raised to Full access, which needs the access level of the sharing user's row in the Share dialog
  and its privilege tree. The words of the recipient's line in the Sharing pane ("has special
  permissions") are claimed on the pane as a whole.

  The recipient of the share is the library's sharing user. Names are letters and digits only (the
  Dashboards search misses "-" and "_") and carry the run's time. The project, under both names, is
  deleted with its table, view and grant when the feature starts and ends. It is serial: the
  Dashboards search is shared with every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeFiles{time}" is on the server
    And no project named "BDDLifeFilesRenamed{time}" is on the server

  Scenario: The file is saved as a project with Data sync
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDLifeFiles{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeFiles{time}" uploaded' should have been shown
    And 1 project named "BDDLifeFiles{time}" should be on the server
    And the "demog" table of the "BDDLifeFiles{time}" project should be saved with data sync
    And "Share BDDLifeFiles{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeFiles{time}" dialog
    Then the "Share BDDLifeFiles{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: The project is shared to view and use, without notifications
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeFiles{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDLifeFiles{time} gallery card
    Then "Share BDDLifeFiles{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeFiles{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDLifeFiles{time}" dialog
    When user clicks on OK button in "Share BDDLifeFiles{time}" dialog
    Then the "Share BDDLifeFiles{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on BDDLifeFiles{time} gallery card
    Then the context panel should show "BDDLifeFiles{time}"
    And the sharing pane should list the sharing user
    And Sharing pane in context panel should contain text "has special permissions"

  Scenario: The project is renamed from its card and opens under its new name
    When user picks "Rename..." from the context menu of BDDLifeFiles{time} gallery card
    Then Rename project dialog should be visible
    And Name input in Rename project dialog should have value "BDDLifeFiles{time}"
    When user enters "BDDLifeFilesRenamed{time}" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDLifeFilesRenamed{time}" should be on the server
    And 0 projects named "BDDLifeFiles{time}" should be on the server
    When user enters "BDDLifeFilesRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeFilesRenamed{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The owner deletes the project from its card
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeFilesRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDLifeFilesRenamed{time} gallery card
    Then "Are you sure?" dialog should contain text 'Delete project "BDDLifeFilesRenamed{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeFilesRenamed{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDLifeFilesRenamed{time} gallery card should be absent
