@journey @serial @realizes:views.projects @realizes:views.space @realizes:sharing.share-dialog
Feature: A project built from a file in a space, shared, and renamed with its space
  A space is created in the Browse tree and demog.csv is copied into it from Browse > Files > Demo;
  the file opened from the space is saved as a project with Data sync, and the project alone is
  shared with the second account from its Dashboards card. Then the space and the project are
  renamed, and both accounts reopen the project after each rename. Translated from the TestTrack case
  Projects/projects-lifecycle-spaces, which reproduces GROK-18345 (the recipient could not open a
  shared, data-synced project whose table comes from a space).

  The second account is the library's sharing user; it is switched to by a session of its own, not
  by the platform's Logout (which would end the session every worker shares). The page reloads on
  every switch (the switch starts the error and balloon floors afresh), so the Browse panel is
  opened again after each one.

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time; the space's namespace is its name, so the file's path reads <space>:Files/demog.csv. The
  project, under both names, and the space, under both names (with the copied file), are deleted
  when the feature starts and ends. It is @serial: the save uploads a table the second account reads
  back.

  Known failure (reproduced 2026-09-25 on localhost and on dev 1.28.0, by UI gestures only and again
  through project.open(), which never resolves): once the space is renamed, the project no longer
  opens for either account. The table's creation script keeps the old namespace
  (OpenFile("<old space>:Files/demog.csv")); the open ends in a "Data loading error" dialog that
  quotes that line ("Offset is outside the bounds of the DataView") and the Dashboards view stays
  current, while the file itself opens from the renamed space. Each of the four opens after the
  rename is a normal scenario that waits for the open to end either way (the demog view with rows
  re-read by data sync, or the dialog) and claims neither; the @known-failure scenario after it reads
  that settled state, and its first claim is the defect itself: no "Data loading error" dialog quotes
  the old space's OpenFile line. The renames themselves are claimed on the server and pass.
  GROK-18345 itself (the recipient opening the project before any rename) passes.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeSpaceProj{time}" is on the server
    And no project named "BDDLifeSpaceProjRenamed{time}" is on the server
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
    And the file "BDDLifeSpace{time}:Files/demog.csv" should be on the server
    And the file "System:DemoFiles/demog.csv" should be on the server
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
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDLifeSpaceProj{time}" into Name text input in "Save project" dialog
    And user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then "demog" project table in "Save project" dialog should contain text 'OpenFile("BDDLifeSpace{time}:Files/demog.csv")'
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeSpaceProj{time}" uploaded' should have been shown
    And 1 project named "BDDLifeSpaceProj{time}" should be on the server
    And the "demog" table of the "BDDLifeSpaceProj{time}" project should be saved with data sync
    And the creation script of the "demog" table of the "BDDLifeSpaceProj{time}" project on the server should contain 'OpenFile("BDDLifeSpace{time}:Files/demog.csv")'
    And "Share BDDLifeSpaceProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeSpaceProj{time}" dialog
    Then the "Share BDDLifeSpaceProj{time}" dialog should close
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: Only the project is shared with the second account
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
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
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    And user double-clicks on BDDLifeSpaceProj{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current

  Scenario: The owner renames the space
    Given user signs in as themselves again
    And the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user picks "Rename..." from the context menu of BDDLifeSpace{time} tree node inside browse tree
    Then Rename project dialog should be visible
    And Name input in Rename project dialog should have value "BDDLifeSpace{time}"
    When user enters "BDDLifeSpaceRen{time}" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then the "Rename project" dialog should close
    And 1 space named "BDDLifeSpaceRen{time}" should be on the server
    And 0 spaces named "BDDLifeSpace{time}" should be on the server

  Scenario: The owner opens the project after the space was renamed
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    And user double-clicks on BDDLifeSpaceProj{time} gallery card
    Then the open should end in the "demog" view or the "Data loading error" dialog

  @known-failure @realizes:GROK-21025
  Scenario: After the space rename, the owner's open does not fail on the old space name
    Then no "Data loading error" dialog should show 'OpenFile("BDDLifeSpace{time}:Files/demog.csv")'
    And the "demog" view should be current
    And the table should have 5850 rows

  Scenario: The second account opens the project after the space was renamed
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    And user double-clicks on BDDLifeSpaceProj{time} gallery card
    Then the open should end in the "demog" view or the "Data loading error" dialog

  @known-failure @realizes:GROK-21025
  Scenario: After the space rename, the second account's open does not fail on the old space name
    Then no "Data loading error" dialog should show 'OpenFile("BDDLifeSpace{time}:Files/demog.csv")'
    And the "demog" view should be current
    And the table should have 5850 rows

  Scenario: The owner renames the project
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProj{time}" into gallery search
    And user picks "Rename..." from the context menu of BDDLifeSpaceProj{time} gallery card
    Then Name input in Rename project dialog should have value "BDDLifeSpaceProj{time}"
    When user enters "BDDLifeSpaceProjRenamed{time}" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDLifeSpaceProjRenamed{time}" should be on the server
    And 0 projects named "BDDLifeSpaceProj{time}" should be on the server
    When user enters "BDDLifeSpaceProjRenamed{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    Then BDDLifeSpaceProjRenamed{time} gallery card should be visible

  Scenario: The owner opens the renamed project
    When user double-clicks on BDDLifeSpaceProjRenamed{time} gallery card
    Then the open should end in the "demog" view or the "Data loading error" dialog

  @known-failure @realizes:GROK-21025
  Scenario: The renamed project's open does not fail on the old space name for the owner
    Then no "Data loading error" dialog should show 'OpenFile("BDDLifeSpace{time}:Files/demog.csv")'
    And the "demog" view should be current
    And the table should have 5850 rows

  Scenario: The second account opens the renamed project
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProjRenamed{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    And user double-clicks on BDDLifeSpaceProjRenamed{time} gallery card
    Then the open should end in the "demog" view or the "Data loading error" dialog

  @known-failure @realizes:GROK-21025
  Scenario: The renamed project's open does not fail on the old space name for the second account
    Then no "Data loading error" dialog should show 'OpenFile("BDDLifeSpace{time}:Files/demog.csv")'
    And the "demog" view should be current
    And the table should have 5850 rows

  Scenario: The owner deletes the project and the space
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeSpaceProjRenamed{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDLifeSpaceProjRenamed{time} gallery card
    Then "Are you sure?" dialog should contain text 'Delete project "BDDLifeSpaceProjRenamed{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeSpaceProjRenamed{time}" should be on the server
    Given Spaces tree node inside browse tree is expanded
    When user picks "Delete Space" from the context menu of BDDLifeSpaceRen{time} tree node inside browse tree
    Then "Are you sure?" dialog should contain text 'Delete space "BDDLifeSpaceRen{time}"?'
    And "Are you sure?" dialog should contain text "This will delete space and its related data"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDLifeSpaceRen{time}" should be on the server
