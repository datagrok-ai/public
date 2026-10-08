@journey @serial @realizes:views.projects @realizes:views.space @realizes:sharing.share-dialog
Feature: A project built from a file in a space, shared, and reopened after the space is renamed
  A space is created in the Browse tree and demog.csv is copied into it from Browse > Files > Demo;
  the file opened from the space is saved as a project with Data sync (its creation script reads
  the space's file), and the project alone is shared with the second account from its Dashboards
  card. Then the space is renamed and the owner reopens the project. Translated from the TestTrack
  case Projects/projects-lifecycle-spaces.

  The second account (the library's sharing user) signs in on the feature's page for the
  recipient's opens of the project (GROK-18345), before and after the space rename.

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time; the space's namespace is its name, so the file's path reads <space>:Files/demog.csv. The
  project and the space, under both its names (with the copied file), are deleted when the feature
  starts and ends. It is serial: the Dashboards search is shared with every feature that saves a
  project.

  Known failure GROK-21025: once the space is renamed, the project no longer opens (the table's
  creation script keeps the old namespace, OpenFile("<old space>:Files/demog.csv")). The ticket and
  the runs of 1.28.0 in September saw the open end in a "Data loading error" dialog quoting that
  line; on localhost 1.28.0 f98337ae1b (2026-10-05) no dialog comes — "Opening project" stays in
  the task bar for over two minutes and the Dashboards view stays current — so the defect has no
  readable signature to claim, and each known failure (the owner's and the second account's open)
  claims only the demog view with its rows re-read by data sync. The scenario after the owner's
  signs the second account in, and the one after the second account's signs the owner back in.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeSpaceProj{time}" is on the server
    And no space named "BDDLifeSpace{time}, BDDLifeSpaceRen{time}" is on the server

  Scenario: A space is created
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDLifeSpace{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDLifeSpace{time}" should be on the server
    Given Spaces tree node inside browse tree is expanded
    Then BDDLifeSpace{time} tree node inside browse tree should be visible

  Scenario: The demo file is copied into the space
    Given Files tree node inside browse tree is expanded
    When user clicks on "Files > Demo" tree node inside browse tree
    Then the "Demo" view should be current
    When user drags demog.csv link in gallery to BDDLifeSpace{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should have value "Link"
    And choice input in Move entity dialog should offer "Move, Link, Copy"
    When user selects "Copy" in Move entity dialog
    Then choice input in Move entity dialog should have value "Copy"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDDLifeSpace{time} tree node inside browse tree
    Then the "BDDLifeSpace{time}" view should be current
    And demog.csv link in gallery should be visible
    When user clicks on "Files > Demo" tree node inside browse tree
    Then the "Demo" view should be current
    And demog.csv link in gallery should be visible

  Scenario: The file opens from the space
    When user double-clicks on BDDLifeSpace{time} tree node inside browse tree
    Then the "BDDLifeSpace{time}" view should be current
    When user double-clicks on demog.csv link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the page address should contain "/file/BDDLifeSpace{time}.Files/demog.csv"

  Scenario: The file from the space is saved as a project with Data sync
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDLifeSpaceProj{time}" into Name text input in "Save project" dialog
    And user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then "demog" project table in "Save project" dialog should contain text 'OpenFile("BDDLifeSpace{time}:Files/demog.csv")'
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeSpaceProj{time}" uploaded' should have been shown
    And 1 project named "BDDLifeSpaceProj{time}" should be on the server
    And the "demog" table of the "BDDLifeSpaceProj{time}" project should be saved with data sync
    And "Share BDDLifeSpaceProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeSpaceProj{time}" dialog
    Then the "Share BDDLifeSpaceProj{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: Only the project is shared with the second account
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDLifeSpaceProj{time} gallery card
    Then "Share BDDLifeSpaceProj{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeSpaceProj{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDLifeSpaceProj{time}" dialog
    And user clicks on OK button in "Share BDDLifeSpaceProj{time}" dialog
    Then the "Share BDDLifeSpaceProj{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on BDDLifeSpaceProj{time} gallery card
    Then the context panel should show "BDDLifeSpaceProj{time}"
    And the sharing pane should list the sharing user

  Scenario: The second account opens the shared project (GROK-18345)
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeSpaceProj{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown
    When user picks "Close All" from the context menu of browse tab
    And user signs in as themselves again
    Then the running account should be signed in

  Scenario: The owner renames the space
    Then the running account should be signed in
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user picks "Rename..." from the context menu of BDDLifeSpace{time} tree node inside browse tree
    Then Rename project dialog should be visible
    And Name input in Rename project dialog should have value "BDDLifeSpace{time}"
    When user enters "BDDLifeSpaceRen{time}" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then the "Rename project" dialog should close
    And 1 space named "BDDLifeSpaceRen{time}" should be on the server
    And 0 spaces named "BDDLifeSpace{time}" should be on the server
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    Then BDDLifeSpaceProj{time} gallery card should be visible

  @known-failure @realizes:GROK-21025
  Scenario: After the space rename, the owner opens the project with its rows
    When user double-clicks on BDDLifeSpaceProj{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync

  Scenario: The second account signs in after the space rename
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then BDDLifeSpaceProj{time} gallery card should be visible

  @known-failure @realizes:GROK-21025
  Scenario: After the space rename, the second account opens the project with its rows
    When user double-clicks on BDDLifeSpaceProj{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync

  Scenario: The owner deletes the project and the space
    When user signs in as themselves again
    Then the running account should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDLifeSpaceProj{time} gallery card
    Then "Are you sure?" dialog should contain text 'Delete project "BDDLifeSpaceProj{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeSpaceProj{time}" should be on the server
    Given Spaces tree node inside browse tree is expanded
    When user picks "Delete Space" from the context menu of BDDLifeSpaceRen{time} tree node inside browse tree
    Then "Are you sure?" dialog should contain text 'Delete space "BDDLifeSpaceRen{time}"?'
    And "Are you sure?" dialog should contain text "This will delete space and its related data"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDLifeSpaceRen{time}" should be on the server
