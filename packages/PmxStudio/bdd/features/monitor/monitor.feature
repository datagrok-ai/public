Feature: Monitor a running study
  Whoever follows a study sees who needs attention, and reveals the groups only when someone must decide;
  revealing them goes in the study's log.

  Background:
    Given user is logged in
    And user opens PMX Studio on Monitor

  @journey
  Scenario: Follow the example study and reveal the groups
    When user clicks on "Follow an example study" button
    Then the study should have finished running
    And study error should be hidden
    And "Needs attention" section should be visible
    And "Measurement sheet" section should be visible
    When user clicks on "Decisions" tab
    And user clicks on "Reveal the groups" button
    Then "Log" section should contain text "Groups revealed"
