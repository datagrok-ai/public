@serial
Feature: A request the user stopped stays stopped
  Stop withdraws a request: a later prompt that asks for something else, or carries no request at all, must
  not resume it. The job adds its columns one at a time, so a few seconds in it cannot have added the last
  one, and that column must still be missing after the follow-ups.

  These scenarios reach a Docker container and a paid model, which the library's hard rule keeps out of
  the regular suite: they run only with GROKKY_BDD_LIVE=1 and skip otherwise. @serial: one chat at a time.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens mol1K dataset
    And the AI panel is open on a new chat

  Scenario: A punctuation-only prompt after Stop does not redo the stopped job
    When user sends "Add 8 columns one by one named stop_1 to stop_8, each equal to Q plus its number, one step each" to the AI panel
    And user waits 6 seconds
    And user clicks on AI stop button
    Then the AI turn should have been stopped
    And the table should not have a column "stop_8"
    When user asks the AI panel "How many rows does the table have?"
    Then the table should not have a column "stop_8"
    When user asks the AI panel "... ??? !!!"
    Then the table should not have a column "stop_8"
    And the stopped AI turn should have sent nothing after it was aborted
    And a screenshot "F-007-after-punctuation" of the AI panel is saved
