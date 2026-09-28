@journey @serial @realizes:views.projects @realizes:views.space
Feature: A project moved into a space and on to another
  Two spaces are created in the Browse tree and demog.csv is saved as a project; the project is moved
  into the first space from its Dashboards card, opened from there, moved on to the second space from
  its link in the first, and opened again. Translated from the TestTrack case Projects/complex-move
  (step 10 of Projects/complex-ui, the move to a space; a file share cannot take an entity).

  Every move is claimed three ways: the balloon, the link in the space's own view, and the server
  (the space lists the project among its children and the project's grok name is in the space's
  namespace); leaving the first space is claimed in its view and on the server.

  After each open the claim is the md's: no error balloon. A double click on the project's link in a
  space, when the project is already the current object (the card was right-clicked for the move, or
  the link clicked once), opens the project and also shows the warning 'Project "<name>" is already
  open' (seen three times on localhost, 2026-09-25; from a fresh current object the same double click
  shows none). That warning is a suspected defect, not claimed here either way.

  GROK-21022 is the known cause of this feature's occasional failure (2 runs of 7 on localhost, core
  1.28.0 bc64f40e47, 2026-09-25; 0 of 9 on 2026-09-28): after the project is opened from a space's
  link and Close All, a double click on a space node in the Browse tree sometimes opens nothing — the
  Home view stays current and the console gets "NullError: method not found ... on null" at
  browse_panel_preview.dart 90:56, `shell.previewProjects.removeWhere((p) => p.id !=
  parentProject.id)` with no parent dashboard (a space node has none) while a preview project is still
  listed. A preview open of a project that another open registers first leaves its copy in that list,
  and Close All does not remove it; that takes two opens of one project racing, so no UI sequence
  reproduces it every time and it is not tagged @known-failure. The claims stand: the space view
  should open.

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The project and both spaces are deleted when the feature starts and ends. It is @serial: the
  save uploads a table the reopen reads back.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDMove{time}" is on the server
    And no space named "BDDMoveA{time}, BDDMoveB{time}" is on the server

  Scenario: Two spaces are created
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDMoveA{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDMoveA{time}" should be on the server
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDDMoveB{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDMoveB{time}" should be on the server
    Given Spaces tree node inside browse tree is expanded
    Then BDDMoveA{time} tree node inside browse tree should be visible
    And BDDMoveB{time} tree node inside browse tree should be visible

  Scenario: The file is saved as a project
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDMove{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDMove{time}" uploaded' should have been shown
    And 1 project named "BDDMove{time}" should be on the server
    And "Share BDDMove{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDMove{time}" dialog
    Then the "Share BDDMove{time}" dialog should close
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no errors should have been logged

  Scenario: The project is moved into the first space from its card
    Given the browse panel is open
    And the "BDDMoveA{time}" space should not hold the project "BDDMove{time}" on the server
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDMove{time}" into gallery search
    And user picks "Move to Space..." from the context menu of BDDMove{time} gallery card
    Then Move to space dialog should be visible
    When user selects "BDDMoveA{time}" in Space input in Move to space dialog
    And user clicks on OK button in Move to space dialog
    Then the "Move to space" dialog should close
    And an info balloon containing "Moved BDDMove{time} to BDDMoveA{time}" should have been shown
    And the "BDDMoveA{time}" space should hold the project "BDDMove{time}" on the server
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current
    And BDDMove{time} link in space gallery should be visible

  Scenario: The project opens from the first space
    When user double-clicks on BDDMove{time} link in space gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And no error balloon should have been shown
    And no errors should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current

  Scenario: The project is moved on to the second space from the first
    Given the browse panel is open
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current
    When user picks "Move to Space..." from the context menu of BDDMove{time} link in space gallery
    Then Move to space dialog should be visible
    When user selects "BDDMoveB{time}" in Space input in Move to space dialog
    And user clicks on OK button in Move to space dialog
    Then the "Move to space" dialog should close
    And an info balloon containing "Moved BDDMove{time} to BDDMoveB{time}" should have been shown
    And the "BDDMoveB{time}" space should hold the project "BDDMove{time}" on the server
    And the "BDDMoveA{time}" space should not hold the project "BDDMove{time}" on the server
    When user double-clicks on BDDMoveB{time} tree node inside browse tree
    Then the "BDDMoveB{time}" view should be current
    And BDDMove{time} link in space gallery should be visible
    When user double-clicks on BDDMoveA{time} tree node inside browse tree
    Then the "BDDMoveA{time}" view should be current
    And BDDMove{time} link in space gallery should be absent

  Scenario: The project opens from the second space
    When user double-clicks on BDDMoveB{time} tree node inside browse tree
    Then the "BDDMoveB{time}" view should be current
    When user double-clicks on BDDMove{time} link in space gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And no error balloon should have been shown
    And no errors should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current

  Scenario: The project and both spaces are deleted
    Given the browse panel is open
    When user double-clicks on BDDMoveB{time} tree node inside browse tree
    Then the "BDDMoveB{time}" view should be current
    When user picks "Delete Project" from the context menu of BDDMove{time} link in space gallery
    Then "Are you sure?" dialog should contain text 'Delete project "BDDMove{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDMove{time}" should be on the server
    When user picks "Delete Space" from the context menu of BDDMoveA{time} tree node inside browse tree
    Then "Are you sure?" dialog should contain text 'Delete space "BDDMoveA{time}"?'
    And "Are you sure?" dialog should contain text "This will delete space and its related data"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDMoveA{time}" should be on the server
    When user picks "Delete Space" from the context menu of BDDMoveB{time} tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDMoveB{time}" should be on the server
