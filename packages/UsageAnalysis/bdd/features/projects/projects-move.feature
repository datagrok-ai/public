@journey @serial @realizes:views.projects @realizes:views.space @realizes:sharing.share-dialog
Feature: A project moved between spaces, and access through a shared space
  Two root spaces and a child of the first are created in the Browse tree; demog.csv from Browse >
  Files > Demo is saved as a project with Data sync and moved into the first space from its
  Dashboards card (Move to Space...), then dragged into the child space in the Browse tree, then
  moved on to the second space from the child space's view (Move to Space...); it opens with its
  data after every move. The first space, not the project, is shared with the second account, and
  the project's Share dialog shows the access as inherited from the space. Translated from the
  TestTrack case Projects/complex-move.

  The second account (the library's sharing user) signs in on the feature's page and opens the project
  from the shared space and from its child space, and no longer finds it after the move to the
  unshared space. A second project, iris saved as BDDMoveStay, is moved into the child space (Move to
  Space... offers root spaces only, so it goes into the first space and is dragged on) and stays
  there: each "no longer sees" claim follows a claim that this project is listed in the same gallery,
  so an absence cannot be read off a gallery that has not filled. Kept without two claims, restored
  with the phrases (see the request document): that the sharing user is listed under "Inherited from"
  in the project's Share dialog (the dialog's grant rows have no element), and that the recipient's
  Save dialog has Save original project disabled (one option of the dialog's radio group has no
  reading).

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The two projects and the three spaces (the child is removed with its root, and is swept by name
  too) are deleted when the feature starts and ends; the space's grant goes with the space. It is
  serial: the Dashboards search is shared with every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDMoveProj{time}" is on the server
    And no project named "BDDMoveStay{time}" is on the server
    And no space named "BDDMoveA{time}, BDDMoveAChild{time}, BDDMoveB{time}" is on the server

  Scenario: Two spaces and a child of the first are created
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDMoveA{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDMoveB{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    Given Spaces tree node inside browse tree is expanded
    Then BDDMoveA{time} tree node inside browse tree should be visible
    And BDDMoveB{time} tree node inside browse tree should be visible
    When user picks "Create Child Space..." from the context menu of BDDMoveA{time} tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDMoveAChild{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And BDDMoveAChild{time} tree node inside browse tree should be visible
    And 1 space named "BDDMoveA{time}" should be on the server
    And 1 space named "BDDMoveB{time}" should be on the server

  Scenario: A file is saved as a project with Data sync
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDMoveProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDMoveProj{time}" uploaded' should have been shown
    And "Share BDDMoveProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDMoveProj{time}" dialog
    Then the "Share BDDMoveProj{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: The project is moved into the first space from its card and opens there
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDMoveProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Move to Space..." from the context menu of BDDMoveProj{time} gallery card
    Then "Move to space" dialog should be visible
    When user selects "BDDMoveA{time}" in Space input in "Move to space" dialog
    And user clicks on OK button in "Move to space" dialog
    Then the "Move to space" dialog should close
    Given Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current
    And "BDDMoveProj{time}" link in gallery should be visible
    When user double-clicks on "BDDMoveProj{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And no error or warning balloon should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current

  Scenario: The space is shared, and the project's access is inherited from it
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of BDDMoveA{time} tree node inside browse tree
    Then "Share BDDMoveA{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDMoveA{time}" dialog
    And user clicks on OK button in "Share BDDMoveA{time}" dialog
    Then the "Share BDDMoveA{time}" dialog should close
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current
    When user picks "Share..." from the context menu of "BDDMoveProj{time}" link in gallery
    Then "Share BDDMoveProj{time}" dialog should be visible
    And "Share BDDMoveProj{time}" dialog should contain text "Inherited from"
    And "Share BDDMoveProj{time}" dialog should contain text "BDDMoveA{time}"
    When user clicks on CANCEL button in "Share BDDMoveProj{time}" dialog
    Then the "Share BDDMoveProj{time}" dialog should close

  Scenario: The second account opens the project from the shared space
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current
    And "BDDMoveProj{time}" link in gallery should be visible
    When user double-clicks on "BDDMoveProj{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And "Data loading error" dialog should be absent
    When user picks "Close All" from the context menu of browse tab
    And user signs in as themselves again
    Then the running account should be signed in
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current

  Scenario: The project is dragged into the child space and opens there
    Then the running account should be signed in
    When user drags "BDDMoveProj{time}" link in gallery to BDDMoveAChild{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDDMoveAChild{time} tree node inside browse tree
    Then the "BDDMoveAChild{time}" view should be current
    And "BDDMoveProj{time}" link in gallery should be visible
    When user double-clicks on "BDDMoveProj{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And no error or warning balloon should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current

  Scenario: The second account opens the project in the child space
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    And BDDMoveA{time} tree node inside browse tree is expanded
    When user double-clicks on BDDMoveAChild{time} tree node inside browse tree
    Then the "BDDMoveAChild{time}" view should be current
    And "BDDMoveProj{time}" link in gallery should be visible
    When user double-clicks on "BDDMoveProj{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And "Data loading error" dialog should be absent
    When user picks "Close All" from the context menu of browse tab
    And user signs in as themselves again
    Then the running account should be signed in

  Scenario: The project is moved on to the second space and opens there
    Then the running account should be signed in
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDMoveAChild{time} tree node inside browse tree
    Then the "BDDMoveAChild{time}" view should be current
    When user picks "Move to Space..." from the context menu of "BDDMoveProj{time}" link in gallery
    Then "Move to space" dialog should be visible
    When user selects "BDDMoveB{time}" in Space input in "Move to space" dialog
    And user clicks on OK button in "Move to space" dialog
    Then the "Move to space" dialog should close
    When user double-clicks on BDDMoveB{time} tree node inside browse tree
    Then the "BDDMoveB{time}" view should be current
    And "BDDMoveProj{time}" link in gallery should be visible
    When user double-clicks on BDDMoveAChild{time} tree node inside browse tree
    Then the "BDDMoveAChild{time}" view should be current
    And "BDDMoveProj{time}" link in gallery should be absent
    When user double-clicks on BDDMoveB{time} tree node inside browse tree
    Then the "BDDMoveB{time}" view should be current
    When user double-clicks on "BDDMoveProj{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And no error or warning balloon should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current

  Scenario: A second project is saved and moved into the child space
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---iris.csv tree node inside browse tree
    Then the "iris" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDMoveStay{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDMoveStay{time}" uploaded' should have been shown
    When user clicks on CANCEL button in "Share BDDMoveStay{time}" dialog
    Then the "Share BDDMoveStay{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDMoveStay{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Move to Space..." from the context menu of BDDMoveStay{time} gallery card
    Then "Move to space" dialog should be visible
    When user selects "BDDMoveA{time}" in Space input in "Move to space" dialog
    And user clicks on OK button in "Move to space" dialog
    Then the "Move to space" dialog should close
    Given Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current
    Given BDDMoveA{time} tree node inside browse tree is expanded
    When user drags "BDDMoveStay{time}" link in gallery to BDDMoveAChild{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDDMoveAChild{time} tree node inside browse tree
    Then the "BDDMoveAChild{time}" view should be current
    And "BDDMoveStay{time}" link in gallery should be visible

  Scenario: After the move to the unshared space, the second account no longer sees the project
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDMove" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then BDDMoveStay{time} gallery card should be visible
    And BDDMoveProj{time} gallery card should be absent
    Given Spaces tree node inside browse tree is expanded
    Then BDDMoveA{time} tree node inside browse tree should be visible
    And BDDMoveB{time} tree node inside browse tree should be absent
    Given BDDMoveA{time} tree node inside browse tree is expanded
    When user double-clicks on BDDMoveAChild{time} tree node inside browse tree
    Then the "BDDMoveAChild{time}" view should be current
    And "BDDMoveStay{time}" link in gallery should be visible
    And "BDDMoveProj{time}" link in gallery should be absent
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: The project and the spaces are deleted
    Then the running account should be signed in
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDMoveB{time} tree node inside browse tree
    Then the "BDDMoveB{time}" view should be current
    When user picks "Delete Project" from the context menu of "BDDMoveProj{time}" link in gallery
    Then "Are you sure?" dialog should contain text 'Delete project "BDDMoveProj{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDMoveProj{time}" should be on the server
    When user picks "Delete Space" from the context menu of BDDMoveA{time} tree node inside browse tree
    Then "Are you sure?" dialog should contain text 'Delete space "BDDMoveA{time}"?'
    And "Are you sure?" dialog should contain text "This will delete space and its related data"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user picks "Delete Space" from the context menu of BDDMoveB{time} tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDMoveA{time}" should be on the server
    And 0 spaces named "BDDMoveB{time}" should be on the server
    And no errors should have been logged
