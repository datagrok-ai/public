@journey
Feature: Analyze
  The Analyze view element by element on its example study: the start page, every tab, an animal picked, the
  tumor-corrected weights, the texts ready to paste, and the form's fields.

  Background:
    Given user is logged in
    And user opens PMX Studio on Analyze

  Scenario: Start page
    Then analyze start should be visible
    And "Analyze an example study" button should be enabled
    And form run button should be disabled

  Scenario: Example study, Efficacy tab
    When user clicks on "Analyze an example study" button
    Then the study should have finished running
    And study error should be hidden
    And analyze start should be hidden
    And stale note should be hidden
    And "Experimental groups" section should contain text "Reference"
    And "Experimental groups" section should contain text "TGI 70%"
    And "Observed tumor growth" section should be visible
    When user clicks on linear axis choice
    Then visible volume axis should have the value "Linear"

  Scenario: Each animal tab
    When user clicks on "Each animal" tab
    Then "Each animal against its untreated twin" section should be visible
    And animal list should contain text "respond"
    And twin effect line should contain text "than its untreated twin"
    When user remembers the twin line
    And user clicks on the last animal row
    Then the twin line should differ from remembered

  Scenario: Tolerability tab
    When user clicks on "Tolerability" tab
    Then "Weight loss" section should be visible
    And "Body weight change from day 0" section should contain text "as weighed"
    When user checks "Tumor-corrected" input
    Then "Body weight change from day 0" section should contain text "without the tumor"
    When user unchecks "Tumor-corrected" input
    Then "Body weight change from day 0" section should contain text "as weighed"

  Scenario: Plan vs observed tab
    When user clicks on "Plan vs observed" tab
    Then "Effect of each treatment" section should be visible
    And "Control growth against the prediction" section should be visible
    And "End of the study" section should contain text "median control"
    And "Contribute to the atlas" section should contain text "Coming soon"

  Scenario: Statistics tab
    When user clicks on "Statistics" tab
    Then "Results summary" section should contain text "evaluated on day"
    And "Growth rates" section should be visible
    And "Control check" section should contain text "control measurements"
    And "Analysis protocol" section should contain text "Primary comparison"
    When user clicks on first copy button
    Then an info balloon containing "Results summary copied" should have been shown
    And the clipboard should contain text "evaluated on day"

  Scenario: The form's fields mark the analysis out of date and Run brings them in
    When user enters "1000" into "Endpoint volume" input
    Then stale note should be visible
    When user remembers the "Evaluation day" figure
    And user clicks on form run button
    Then the study should have finished running
    And stale note should be hidden
    And the "Evaluation day" figure should be lower than remembered
    When user enters "100" into "Endpoint volume" input
    Then form run button should be disabled
    When user clears "Endpoint volume" input
    Then form run button should be disabled
    When user enters "1500" into "Endpoint volume" input
    And user selects "CT26" in "Cell line" input
    And user clicks on form run button
    Then the study should have finished running
    And study error should be hidden
    When user selects "MC38" in "Cell line" input
    And user clicks on form run button
    Then the study should have finished running

  Scenario: Advanced fields
    When user clicks on advanced header
    Then "Evaluation day" input should be visible
    And "Significance (α)" input should be visible
    And "Humane limit, weight loss" input should be visible
    And "Humane limit, tumor" input should be visible
    And "Tumor density" input should be visible
    When user enters "9" into "Evaluation day" input
    And user clicks on form run button
    Then the study should have finished running
    And the "Evaluation day" figure should read "Day 9"
    When user clears "Evaluation day" input
    And user enters "1.2" into "Tumor density" input
    And user clicks on form run button
    Then the study should have finished running
    And study error should be hidden

  Scenario: An evaluation day that cannot be analyzed says why
    When user enters "5" into "Evaluation day" input
    And user clicks on form run button
    Then the study should have finished running
    And study error should contain text "choose a measurement day"
    When user enters "0" into "Evaluation day" input
    And user clicks on form run button
    Then the study should have finished running
    And study error should contain text "choose a measurement day"
    When user clears "Evaluation day" input
    And user clicks on form run button
    Then the study should have finished running
    And study error should be hidden
