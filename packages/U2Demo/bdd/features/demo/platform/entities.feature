@demo @realizes:u2.dg.entities
Feature: Entity pickers and chips
  The pickers are inline typeaheads without labels, so they go by the placeholder the user sees.

  Scenario: Picking a user and a group renders them as chips
    Given user opens the "Entities" demo page
    Then Admin chip should be visible
    When user types "adm" into "Search users…" typeahead
    And user selects "Admin" in "Search users…" typeahead
    Then value of user readout should have text "Admin"
    When user types "All users" into "Search groups…" typeahead
    And user selects "All users" in "Search groups…" typeahead
    Then value of group readout should have text "All users"
    And "All users" chip should be visible
