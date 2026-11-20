@serial
Feature: A prompt sent while the AI panel answers is not silently dropped
  Enter in the prompt box while a turn runs does not start a second turn; the panel must tell the user
  to wait for the response or press Stop instead of doing nothing visible.

  These scenarios reach a Docker container and a paid model: they run only with GROKKY_BDD_LIVE=1.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens mol1K dataset
    And the AI panel is open on a new chat

  Scenario: Enter on a second prompt during a turn tells the user to wait
    Given the notification balloons are closed
    When user sends "In your terminal, run `sleep 61 && echo finished` and tell me what it printed." to the AI panel
    And user presses Enter on "Also tell me how many rows the table has" in the AI panel while it answers
    Then an info balloon containing "Wait for the response or press Stop" should have been shown
    When user clicks on AI stop button
    Then the AI turn should have been stopped
