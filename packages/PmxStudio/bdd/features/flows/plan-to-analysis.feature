@journey
Feature: From the plan to the analysis
  A study planned and saved, then analyzed against that plan from Analyze's list of saved plans.

  Scenario: The saved plan is what the analysis compares with
    Given user is logged in
    And no plan named "BDD flow" is saved
    And user opens PMX Studio
    Then the study should have finished running
    When user clicks on "Save plan" button
    And user enters "BDD flow" into "Name" input in "Save plan" dialog
    And user clicks on "Save" button in "Save plan" dialog
    Then a plan named "BDD flow" should be saved
    When user clicks on "3 Analyze" button
    Then analyze start should be visible
    When user clicks on "Analyze an example study" button
    Then the study should have finished running
    When user picks the saved plan "BDD flow"
    Then stale note should be visible
    When user clicks on form run button
    Then the study should have finished running
    And study error should be hidden
    And stale note should be hidden
    When user clicks on "Plan vs observed" tab
    Then results pane should contain text "plan: BDD flow"
    And "Effect of each treatment" section should be visible
