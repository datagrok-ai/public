@serial
Feature: The AI panel's chat survives a stop, a bash command, a reload and a lost agent process
  A chat is one runtime session: a turn stopped mid-run must end everything it started, a !command runs in
  bash mode for that prompt alone, and the chat must go on after the browser reloads and after the runtime
  loses the chat's agent process. The follow-ups only make sense with the earlier turn's answer, and every
  answer is checked against the value computed from the table in the page; the runtime's side is read from
  its WebSocket traffic and, through `docker exec` on this machine, from its container.

  These scenarios reach a Docker container and a paid model, which the library's hard rule keeps out of
  the regular suite: they run only with GROKKY_BDD_LIVE=1 and skip otherwise. @serial: one chat at a time.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens demog-1000 dataset
    And the AI panel is open on a new chat

  Scenario: A long job stopped mid-run leaves nothing running, and the chat goes on
    When user sends "In your terminal, run `sleep 613 && echo finished` (it stands in for a long data job) and tell me what it printed." to the AI panel
    Then the AI runtime container should be running "sleep 613"
    When user clicks on AI stop button
    Then the AI turn should have been stopped
    And the AI panel should show the message "Processing aborted by user"
    And the AI runtime container should not be running "sleep 613"
    When user asks the AI panel "What is the median WEIGHT of the patients in this table, to one decimal? Reply with just the number."
    Then the AI answer should be the median of "WEIGHT"
    And the AI should have used one of the tools "datagrok_exec, terminal"
    And the stopped AI turn should have sent nothing after it was aborted
    And every turn of the chat should have used the same AI session
    And the AI runtime should have kept the chat's history itself
    And a screenshot "session-stop" of the AI panel is saved

  Scenario: A !command shows its raw output, and the next prompt is back in full mode
    When user asks the AI panel "!hostname"
    Then the AI answer should be exactly the host name of the AI runtime container
    And the AI should have used the "terminal" tool
    And the last prompt should have reached the runtime in bash mode
    When user asks the AI panel "What is the median AGE of the patients in this table, to one decimal? Reply with just the number."
    Then the AI answer should be the median of "AGE"
    And the last prompt should have reached the runtime in full mode
    And every turn of the chat should have used the same AI session
    And a screenshot "session-bash" of the AI panel is saved

  Scenario: The chat goes on after the browser reloads
    When user asks the AI panel "What is the median HEIGHT of the patients in this table, in cm, to one decimal? Reply with just the number."
    Then the AI answer should be the median of "HEIGHT"
    When user reloads the page
    And user opens demog-1000 dataset
    And user reopens the AI chat from the chat history
    Then the AI panel should show the earlier turns of the chat
    When user asks the AI panel "What is that in inches, to one decimal? Reply with just the number."
    Then the AI answer should be the median of "HEIGHT" divided by 2.54
    And every turn of the chat should have used the same AI session
    And the AI runtime should have kept the chat's history itself
    And a screenshot "session-reload" of the AI panel is saved

  Scenario: The chat goes on after the runtime loses the chat's agent process
    When user asks the AI panel "What is the median HEIGHT of the patients in this table, in cm, to one decimal? Reply with just the number."
    Then the AI answer should be the median of "HEIGHT"
    And the AI runtime container should run the Hermes agent and not the Claude Agent SDK
    When the AI runtime's agent processes are killed
    And user asks the AI panel "What is that in inches, to one decimal? Reply with just the number."
    Then the AI answer should be the median of "HEIGHT" divided by 2.54
    And every turn of the chat should have used the same AI session
    And the AI runtime should have kept the chat's history itself
    And a screenshot "session-agent-killed" of the AI panel is saved
