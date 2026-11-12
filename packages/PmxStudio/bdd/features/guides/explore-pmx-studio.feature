@guide
Feature: Explore PMX Studio
  A guide: the answer to "what does PMX Studio do, and how do I find my way around it?". The Overview lists the
  lab's studies; Plan sizes one from the atlas of historical controls, Monitor follows it while it runs, Analyze reads
  it once it has run; the header moves between them. The study followed is the example one (a published MC38 table),
  opened with the plan saved in Plan.

  Scenario: Plan a study, follow it and analyze it
    Given the guide may take up to 15 minutes
    And user is logged in
    And simple mode is off
    And the stages show their default layouts
    And no plan named "MC38 guide plan" is saved
    And the shell shows no tooltips
    And the browse panel is open
    # caption: Expand Apps in the Browse panel
    When user expands Apps tree node inside browse tree
    # caption: Open PMX Studio: the Overview lists your studies
    And user clicks on Apps---PMX-Studio tree node inside browse tree
    Then open a study strip should be visible
    # caption: Start a new study in Plan: the tumor model comes first
    When user clicks on plan card and the stage opens
    # caption: Add the treatment arms the study compares with its control
    And user clicks on "Add treatment arm" button and PMX Studio settles
    And user clicks on "Add treatment arm" button and PMX Studio settles
    And user clicks on "Add treatment arm" button and PMX Studio settles
    # caption: The plan: each group's predicted growth, from historical control studies
    Then "Predicted tumor growth" section should be visible
    # caption: Change the endpoint volume: the design follows
    When user enters "1200" into "Endpoint volume" input and PMX Studio settles
    # caption: Scroll to the study calendar
    And user scrolls the page to study calendar
    # caption: Save the plan
    When user clicks on "Save plan" button
    And user enters "MC38 guide plan" into "Name" input
    And user clicks on "Save" button in "Save plan" dialog
    Then a plan named "MC38 guide plan" should be saved
    # caption: Back on the Overview, open a study: its plan, then its measurements so far
    When user clicks on overview step and the stage opens
    And user clicks on open a study strip
    And user clicks on saved plans link
    And user clicks on guide plan menu item
    # caption: The measurements: here an example table, a published MC38 study
    And user uploads "../files/example/mc38-ld4172-anti-pd1.csv" through second file zone
    # caption: Still running: follow it in Monitor
    And user clicks on still running choice
    And user clicks on "Open in Monitor" button in "Open a study" dialog and the stage opens
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
