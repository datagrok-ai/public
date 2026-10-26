@edge-cases
Feature: Table Manager edge cases: no table open
  The Table Manager lists every open table (table-manager.feature, TestTrack General/table-manager.md),
  so with no table open Alt+T docks it with an empty list and no error, and Alt+T closes it again.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And no table view is open
    Then no table should be left in the workspace

  Scenario: With no table open the manager docks with an empty list and no error
    Then "Tables" dock panel should be absent
    When user presses Alt+T
    Then "Tables" dock panel should be visible
    And the "rows" reading of Grid viewer in "Tables" dock panel should be 0
    And no errors should have been logged
    When user presses Alt+T
    Then "Tables" dock panel should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
