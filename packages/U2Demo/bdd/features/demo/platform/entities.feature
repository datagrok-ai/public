@demo @realizes:u2.dg.entities
Feature: Entity pickers and chips
  The pickers are inline typeaheads without labels, so the package names them on its context.

  Scenario: Picking a user and a group renders them as chips
    Given user opens the "Entities" demo page
    Then Admin chip should be visible
    When user types "adm" into user picker
    And user selects "Admin" in user picker
    Then value of user readout should have text "Admin"
    When user types "All users" into group picker
    And user selects "All users" in group picker
    Then value of group readout should have text "All users"
    And "All users" chip should be visible
