@journey @serial @realizes:views.projects @realizes:sharing.share-dialog
Feature: A file-based project shared at two access levels and renamed
  demog.csv opened from Browse > Files > Demo is saved as a project, shared from its Dashboards card
  with the second account to view and use, then raised to full access, and renamed. The second
  account opens it at every stage from its own Dashboards gallery: with View and use its Save dialog
  offers only a copy, with Full access it saves the original. Translated from the TestTrack case
  Projects/projects-lifecycle-files.

  The second account is the library's sharing user; it is switched to by a session of its own, not
  by the platform's Logout (which would end the session every worker shares). The page reloads on
  every switch (the switch starts the error and balloon floors afresh), so the Browse panel is
  opened again after each one. The md lists "Write access"
  among the privileges of the tree; a project's tree has none (Full access, View and use, Edit,
  Delete, Share), so the five that are there are claimed.

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The project, under both names, is deleted when the feature starts and ends, with its table,
  view and grants. It is @serial: the saves upload a table the second account reads back.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeFiles{time}" is on the server
    And no project named "BDDLifeFilesRenamed{time}" is on the server

  Scenario: The file is saved as a project
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDLifeFiles{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeFiles{time}" uploaded' should have been shown
    And 1 project named "BDDLifeFiles{time}" should be on the server
    And "Share BDDLifeFiles{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeFiles{time}" dialog
    Then the "Share BDDLifeFiles{time}" dialog should close
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: The owner shares it to view and use
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeFiles{time}" into gallery search
    And user picks "Share..." from the context menu of BDDLifeFiles{time} gallery card
    Then "Share BDDLifeFiles{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeFiles{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDLifeFiles{time}" dialog
    And user clicks on OK button in "Share BDDLifeFiles{time}" dialog
    Then the "Share BDDLifeFiles{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on BDDLifeFiles{time} gallery card
    Then the context panel should show "BDDLifeFiles{time}"
    And the sharing pane should list the sharing user
    And the sharing pane should show the sharing user as "has special permissions"

  Scenario: With View and use the second account opens it and can only save a copy
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeFiles{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    And user double-clicks on BDDLifeFiles{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And no error or warning balloon should have been shown
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And "Save original project" radio choice in "Save project" dialog should be disabled
    And "Save a copy" radio choice in "Save project" dialog should be checked
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: The owner raises the share to full access
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeFiles{time}" into gallery search
    And user picks "Share..." from the context menu of BDDLifeFiles{time} gallery card
    Then "Share BDDLifeFiles{time}" dialog should be visible
    And the access level of the sharing user in "Share BDDLifeFiles{time}" dialog should be "View and use"
    When user opens the access level of the sharing user in "Share BDDLifeFiles{time}" dialog
    Then the following elements should be visible:
      | Full-access tree node inside privilege tree  |
      | View-and-use tree node inside privilege tree |
      | Edit tree node inside privilege tree         |
      | Delete tree node inside privilege tree       |
      | Share tree node inside privilege tree        |
    And Full-access tree node inside privilege tree should be unchecked
    When user checks Full-access tree node inside privilege tree
    And user clicks outside the privilege tree
    Then the access level of the sharing user in "Share BDDLifeFiles{time}" dialog should be "Full access"
    When user clicks on OK button in "Share BDDLifeFiles{time}" dialog
    Then the "Share BDDLifeFiles{time}" dialog should close
    And an info balloon containing "Shared" should have been shown

  Scenario: The owner renames the project and it still opens
    When user picks "Rename..." from the context menu of BDDLifeFiles{time} gallery card
    And user enters "BDDLifeFilesRenamed{time}" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDLifeFilesRenamed{time}" should be on the server
    And 0 projects named "BDDLifeFiles{time}" should be on the server
    When user enters "BDDLifeFilesRenamed{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    And user double-clicks on BDDLifeFilesRenamed{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: With Full access the second account opens the renamed project and saves the original
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeFiles" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    Then BDDLifeFilesRenamed{time} gallery card should be visible
    When user double-clicks on BDDLifeFilesRenamed{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user remembers when the "BDDLifeFilesRenamed{time}" project was saved
    And user clicks on Save button
    Then "Save project" dialog should be visible
    And "Save original project" radio choice in "Save project" dialog should be enabled
    And "Save original project" radio choice in "Save project" dialog should be checked
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeFilesRenamed{time}" uploaded' should have been shown
    And 1 project named "BDDLifeFilesRenamed{time}" should be on the server
    And the "BDDLifeFilesRenamed{time}" project should have been saved again since remembered
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: The owner deletes the project from its card
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeFilesRenamed{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDLifeFilesRenamed{time} gallery card
    Then "Are you sure?" dialog should contain text 'Delete project "BDDLifeFilesRenamed{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeFilesRenamed{time}" should be on the server
