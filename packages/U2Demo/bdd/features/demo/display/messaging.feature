@demo @realizes:u2.message-input
Feature: The message input

  Scenario: Composing and sending
    Given user opens the "Message input" demo page
    Then Send button should be disabled
    When user types "Hello from bdd" into message input
    Then Send button should be enabled
    When user presses Control+Enter in message input
    Then value of sent readout should contain text "Hello from bdd"
