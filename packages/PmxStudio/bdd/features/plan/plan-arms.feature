@journey
Feature: Plan's treatment arms
  The arm cards: name, therapy, metric and expected effect each reach the plan, arms come and go up to the limit, and
  the multiple-comparisons box follows the number of arms.

  Background:
    Given user is logged in
    And user opens PMX Studio
    Then the study should have finished running

  Scenario: One arm needs no correction
    Then "Correct for multiple comparisons" input should be disabled
    And "Add treatment arm" button should be enabled

  Scenario: Name, therapy, metric and effect reach the plan
    When user enters "Compound A" into "Name" input
    Then stale note should be visible
    When user selects "Cell therapy" in "Therapy" input
    And user clicks on tc choice
    And user enters "40" into "Expected" input
    And user clicks on form run button
    Then the study should have finished running
    And "Experimental groups" section should contain text "Compound A"
    And "Experimental groups" section should contain text "T/C 40%"
    When user enters "99" into "Expected" input
    Then form run button should be disabled
    When user clears "Expected" input
    Then form run button should be disabled
    When user enters "40" into "Expected" input
    Then form run button should be enabled

  Scenario: Arms are added up to the limit and removed
    When user clicks on "Add treatment arm" button
    Then "Correct for multiple comparisons" input should be enabled
    And form run button should be enabled
    And "Multiple comparisons" input should contain text "Dunnett"
    When user clicks on form run button
    Then the study should have finished running
    And "Experimental groups" section should contain text "Treatment 2"
    When user clicks on "Add treatment arm" button
    And user clicks on "Add treatment arm" button
    And user clicks on "Add treatment arm" button
    And user clicks on "Add treatment arm" button
    Then "Add treatment arm" button should be disabled
    When user clicks on form run button
    Then the study should have finished running
    And "Experimental groups" section should contain text "Treatment 6"

  Scenario: Correction off
    When user unchecks "Correct for multiple comparisons" input
    Then "Multiple comparisons" input should contain text "exploratory"
    When user clicks on form run button
    Then the study should have finished running
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "None: each of the 6 comparisons"
    When user clicks on "Design" tab
    And user checks "Correct for multiple comparisons" input
    And user clicks on form run button
    Then the study should have finished running
