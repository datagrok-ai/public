@serial @spaces @edge-cases @realizes:views.space
Feature: Space edge cases: Escape on its dialogs, and a delete that takes a linked project with it
  Guards for what the Create Space, Rename and delete confirmation dialogs keep when they are dismissed
  with Escape, and for what deleting a space removes. Escape closes a dialog as CANCEL does
  (`Dialog._closeOnEsc` in d4/lib/src/widgets/dialog/dialog.dart closes it as canceled); spaces-create,
  spaces-rename and spaces-delete claim CANCEL. A space created through the dialog reaches the server
  seconds after OK, so the space Escape must not create is claimed absent only after a second space,
  created through the same dialog, has arrived.

  Deleting a space removes the space and the links it holds, not a project that was only linked into
  it (core#976: "only the space and the link are deleted"). Today the linked project is deleted with
  the space; that scenario is a known failure, and the tag goes when core#976 is fixed. What a delete
  should do with an entity the deleting user may not delete on its own is not defined (core#378): see
  MISSING.md.

  The spaces and the project are this feature's own, named with the run's time, made through the API
  (the Create Space dialog is spaces-create's subject) and swept when the feature starts and ends.
  The project's name is letters and digits only: the Dashboards search misses "-" and "_".

  Serial with the other features that drag onto or list the Spaces tree.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}" is on the server
    And Spaces tree node inside browse tree is expanded

  Scenario: Escape on the Create Space dialog creates nothing
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDEdgeEsc{time}" into Name input in Create Space dialog
    And user presses Escape
    Then the "Create Space" dialog should close
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDDEdgeMade{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDEdgeMade{time}" should be on the server
    And 0 spaces named "BDDEdgeEsc{time}" should be on the server
    And no errors should have been logged

  Scenario: Escape on the Rename dialog keeps the old name
    Given a space named "BDDEdgeKeep{time}" is on the server
    And Spaces tree node inside browse tree is expanded
    When user picks "Rename..." from the context menu of BDDEdgeKeep{time} tree node inside browse tree
    Then Name input in Rename project dialog should have value "BDDEdgeKeep{time}"
    When user enters "BDDEdgeKeepNew{time}" into Name input in Rename project dialog
    And user presses Escape
    Then Rename project dialog should be hidden
    And 1 space named "BDDEdgeKeep{time}" should be on the server
    And 0 spaces named "BDDEdgeKeepNew{time}" should be on the server
    And BDDEdgeKeep{time} tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Escape on the delete confirmation keeps the space
    Given a space named "BDDEdgeKeep{time}" is on the server
    And Spaces tree node inside browse tree is expanded
    When user picks "Delete Space" from the context menu of BDDEdgeKeep{time} tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    When user presses Escape
    Then the "Are you sure?" dialog should close
    And 1 space named "BDDEdgeKeep{time}" should be on the server
    And BDDEdgeKeep{time} tree node inside browse tree should be visible
    And no errors should have been logged

  @known-failure:94 @realizes:core#976
  Scenario: Deleting a space keeps a project that was only linked into it
    core#976: the project linked into the space is deleted with the space. Remove the tag when
    core#976 is fixed.
    Given no project named "bddedgelinkproj{time}" is on the server
    And a space named "BDDEdgeLink{time}" is on the server
    And user opens demog-1000 dataset
    When user saves the current view as project "bddedgelinkproj{time}"
    And user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "bddedgelinkproj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then bddedgelinkproj{time} gallery card should be visible
    Given Spaces tree node inside browse tree is expanded
    When user drags bddedgelinkproj{time} gallery card to BDDEdgeLink{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should have value "Link"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDDEdgeLink{time} tree node inside browse tree
    Then the "BDDEdgeLink{time}" view should be current
    And "bddedgelinkproj{time}" link in gallery should be visible
    And 1 project named "bddedgelinkproj{time}" should be on the server
    When user picks "Delete Space" from the context menu of BDDEdgeLink{time} tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDEdgeLink{time}" should be on the server
    And 1 project named "bddedgelinkproj{time}" should be on the server
