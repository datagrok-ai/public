@guide
Feature: PMX Studio
  A teaser: the answer to "what does PMX Studio do for a preclinical efficacy study?". The lab's portfolio on the
  Overview; an MC38 study planned, followed halfway through (an animal flagged) and analyzed, the study followed past
  its last measurement and each animal against its digital twin; moving between the stages from the header.

  Scenario: An efficacy study from its design to its analysis
    Given the guide may take up to 10 minutes
    And user is logged in
    And simple mode is off
    And the stages show their default layouts
    And the shell is in presentation mode without tooltips
    And user opens PMX Studio on its Overview
    And the demo study is open in Monitor and Analyze
    # caption: Every tumor efficacy study of the lab, from its design to its analysis
    Then the lab's studies should be on the Overview
    # caption: Planned, running and finished studies side by side, each one a click away
    When user points out your studies
    # caption: Plan: group sizes, duration and power, predicted from historical control studies
    When user clicks on demo plan row and the stage opens
    # caption: The design starts from the tumor model: its growth and variability come from past control groups
    And user points out tumor model section
    # caption: Every assumption is live: lower the endpoint volume and the design follows
    And user drags the slider of "Endpoint volume" input to 1200 and PMX Studio settles
    # caption: Animals per group, study length and power, recomputed at once
    And user points out plan figures
    # caption: The arms to compare, each with its dose and schedule
    And user points out treatment arms section
    # caption: Monitor: the running study against what its plan predicted
    When user clicks on monitor step and the stage opens
    # caption: Each group against its forecast, who needs attention today, and when the study should end
    Then "Needs attention" section should be visible
    # caption: For every animal, the range its next measurement should fall in
    When user clicks on "Measurements" tab
    # caption: Any animal can be flagged for review, with its reason
    And user picks "Flag this animal…" from the context menu of first sheet row
    # caption: Any animal can be flagged for review, with its reason
    And user enters "Ulcerated tumor: check at the next session" into "Reason" input
    # caption: Flagged: the animal carries its mark, in the study table too
    And user clicks on "OK" button and PMX Studio settles
    # caption: Welfare: body weight and condition of every animal, against the humane endpoints
    When user clicks on "Welfare" tab and PMX Studio settles
    # caption: Decisions: flags and exclusions logged with who, when and why; they travel with the study
    When user clicks on "Decisions" tab and PMX Studio settles
    # caption: Analyze: the finished study, every arm against the control
    When user clicks on analyze header step and the stage opens
    # caption: Each group's growth and its effect against the control, with its uncertainty
    And user points out experimental groups section
    # caption: Time to the endpoint as measured, then projected past the last measurement, marked as a projection
    And user scrolls the page to km chart
    # caption: Time to the endpoint as measured, then projected past the last measurement, marked as a projection
    And user drags the slider of km slider to 40
    # caption: Each animal: its response, from the best to the worst, and its tumor over time
    And user scrolls the page to each animal tab
    # caption: Each animal: its response, from the best to the worst, and its tumor over time
    When user clicks on each animal tab and PMX Studio settles
    # caption: Each animal: its response, from the best to the worst, and its tumor over time
    And user points out waterfall chart
    # caption: Its digital twin: the same tumor, predicted as if untreated
    And user scrolls the page to twin section
    # caption: Its digital twin: the same tumor, predicted as if untreated
    And user clicks on third animal row and PMX Studio settles
    # caption: The gap between the animal and its twin is the effect of the treatment, animal by animal
    And user points out twin effect line
    # caption: Tolerability: body weight loss and deaths in every arm
    And user scrolls the page to tolerability tab
    # caption: Tolerability: body weight loss and deaths in every arm
    When user clicks on "Tolerability" tab and PMX Studio settles
    # caption: Against the plan: what the study found against what its design predicted
    When user clicks on "Against the plan" tab and PMX Studio settles
    # caption: The whole portfolio in one place: planned, running and finished
    When user clicks on overview step and the stage opens
