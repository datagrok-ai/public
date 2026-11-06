@serial @spaces @realizes:views.space
Feature: The known failures of spaces
  Temporary: the known failures of spaces while their tickets are open. Each scenario is expected to
  fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario back into
  the feature its description names (spaces-delete.feature) as an ordinary @edge-cases scenario.

  A space deleted is expected to keep a project that was only linked into it (core#976). The space
  and the project are the scenario's own, named with the
  run's time (the project's letters and digits only: the Dashboards search misses "-" and "_"), and
  swept with what they hold when the feature starts and ends.

  Serial with spaces-create, as spaces-delete is: spaces-create making root spaces on the other worker
  while this feature's tree listed Spaces left a fixture space out of the tree.

  Background:
    Given user is logged in
    And the browse panel is open

  @edge-cases @known-failure:48 @realizes:core#976
  Scenario: Deleting a space keeps a project that was only linked into it
    core#976: the project linked into the space is deleted with the space. Remove the tag when
    core#976 is fixed.
    Back into spaces-delete.feature as an ordinary @edge-cases scenario once fixed.
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
