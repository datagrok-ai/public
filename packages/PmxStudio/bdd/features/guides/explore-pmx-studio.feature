@guide
Feature: Explore PMX Studio
  A guide: the answer to "what does PMX Studio do, and how do I find my way around it?". Plan sizes an
  efficacy study from the atlas of historical controls; Analyze reads the study once it has run.

  Scenario: Plan a study and analyze an example
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    # caption: Expand Apps in the Browse panel
    When user expands Apps tree node inside browse tree
    # caption: Open PMX Studio
    And user clicks on Apps---PMX-Studio tree node inside browse tree
    Then the study should have finished running
    And "Experimental groups" section should contain text "Control"
    # caption: Choose the CT26 cell line
    When user selects "CT26" in "Cell line" input
    # caption: Run the plan for the new cell line
    And user clicks on form run button
    Then the study should have finished running
    # caption: See the study calendar
    When user clicks on "Calendar" tab
    Then "Study calendar" section should be visible
    # caption: See how power changes with the design
    When user clicks on "Trade-offs" tab
    Then "Power by group size" section should be visible
    # caption: Read the protocol, ready to paste
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "CT26"
    # caption: Go to the Analyze stage
    When user clicks on "3 Analyze" button
    # caption: Analyze an example study
    And user clicks on "Analyze an example study" button
    Then the study should have finished running
    And "Observed tumor growth" section should be visible
    # caption: Compare each animal with its untreated twin
    When user clicks on "Each animal" tab
    Then "Each animal against its untreated twin" section should be visible
    # caption: Check tolerability
    When user clicks on "Tolerability" tab
    Then "Weight loss" section should be visible
    # caption: Compare the study with its plan
    When user clicks on "Plan vs observed" tab
    Then "Effect of each treatment" section should be visible
    # caption: See the statistics
    When user clicks on "Statistics" tab
    Then "Results summary" section should be visible
