@demo @realizes:u2.object-form
Feature: The object form
  propertyForm() over a real DG.Group: SAVE persists it through dapi and Delete removes it again,
  so the scenario leaves nothing behind on the stand.

  Scenario: Saving and deleting a group through the generated form
    Given user opens the "Object form" demo page
    Then Name input in object form should not be empty
    And value of result readout should have text "(not saved)"
    When user checks Personal checkbox in object form
    And user clicks on SAVE button in object form
    Then value of result readout should contain text "Saved:"
    When user clicks on Delete button in object form
    Then value of result readout should have text "Deleted"
