@demo @realizes:u2.typeahead @realizes:u2.dg.entities
Feature: Space picker
  The space picked in the Space typeahead shows in the preview as the platform's gallery card,
  named by its title, and an entity chip.

  Scenario: A space picked by its name shows in the preview
    Given user is logged in
    And a space named "BDD-U2Demo-Space" is on the server
    And user opens the "Spaces" demo page
    When user types "BDD-U2Demo" into Space typeahead
    And user selects "BDD-U2Demo-Space" in Space typeahead
    Then "BDD-U2Demo-Space" entity card should be visible
    And "BDD-U2Demo-Space" chip should be visible
