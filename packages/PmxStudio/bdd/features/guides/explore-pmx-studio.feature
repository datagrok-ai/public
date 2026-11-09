@guide
Feature: Explore PMX Studio
  A guide: the answer to "what does PMX Studio do, and how do I find my way around it?". The Overview lists the
  lab's studies; Plan sizes one from the atlas of historical controls, Monitor follows it while it runs, Analyze reads
  it once it has run; the header moves between them.

  Scenario: Plan a study, follow it and analyze it
    Given user is logged in
    And simple mode is off
    And the stages show their default layouts
    And the shell shows no tooltips
    And the browse panel is open
    # caption: Expand Apps in the Browse panel
    When user expands Apps tree node inside browse tree
    # caption: Open PMX Studio: the Overview lists your studies
    And user clicks on Apps---PMX-Studio tree node inside browse tree
    Then the studies should be listed
    # caption: Open a planned study
    When user clicks on demo plan row and the stage opens
    # caption: The plan: each group's predicted growth, from historical control studies
    Then "Predicted tumor growth" section should be visible
    # caption: Change the endpoint volume: the design follows
    When user enters "1200" into "Endpoint volume" input and PMX Studio settles
    # caption: Scroll to the study calendar
    And user scrolls the page to study calendar
    # caption: Go to Monitor from the header
    When user clicks on monitor step and the stage opens
    # caption: Follow an example study
    And user clicks on "Follow an example study" button and PMX Studio settles
    # caption: The study today: each group against its forecast, and who needs attention
    Then "Needs attention" section should be visible
    # caption: See each animal's measurements and the range expected next
    When user clicks on "Measurements" tab
    # caption: Check body weight and condition
    When user clicks on "Welfare" tab and PMX Studio settles
    # caption: Go to Analyze from the header: the study you follow, analyzed as of today
    When user clicks on analyze header step and the stage opens
    # caption: Each treated group against the control
    Then "Experimental groups" section should be visible
    # caption: Compare each animal with its untreated twin
    When user clicks on each animal tab and PMX Studio settles
    # caption: The gap between each animal and its twin is the effect of the treatment
    Then "Each animal against its untreated twin" section should be visible
    # caption: Check tolerability
    When user clicks on "Tolerability" tab and PMX Studio settles
    # caption: Compare the study with its plan
    When user clicks on "Against the plan" tab and PMX Studio settles
    # caption: What the study found against what its plan predicted
    Then "Plan vs observed" section should be visible
    # caption: Back to the Overview
    When user clicks on overview step and the stage opens
