@guide
Feature: PMX Studio
  A teaser: the answer to "what does PMX Studio do for a preclinical efficacy study?". The lab's portfolio on the
  Overview; an MC38 study planned, then the example study (a published MC38 table) followed with that plan (an animal flagged) and analyzed, the study followed past
  its last measurement and each animal against its digital twin; moving between the stages from the header.

  Scenario: An efficacy study from its design to its analysis
    Given the guide may take up to 15 minutes
    And user is logged in
    And simple mode is off
    And the stages show their default layouts
    And the shell is in presentation mode without tooltips
    And no plan named "MC38 guide plan" is saved
    And user opens PMX Studio on its Overview
    # caption: Every efficacy study of the lab in one place: planned, running and finished
    When user points out your studies
    # caption: Plan: the design starts from the tumor model and its historical control studies
    When user clicks on plan card and the stage opens
    # caption: Plan: the design starts from the tumor model and its historical control studies
    And user points out tumor model section
    # caption: Three treatment arms, each against the control
    And user clicks on "Add treatment arm" button and PMX Studio settles
    # caption: Three treatment arms, each against the control
    And user clicks on "Add treatment arm" button and PMX Studio settles
    # caption: Three treatment arms, each against the control
    And user clicks on "Add treatment arm" button and PMX Studio settles
    # caption: Every assumption is live: lower the endpoint and group sizes, duration and power follow
    And user drags the slider of "Endpoint volume" input to 1200 and PMX Studio settles
    # caption: Every assumption is live: lower the endpoint and group sizes, duration and power follow
    And user points out plan figures
    # caption: Every assumption is live: lower the endpoint and group sizes, duration and power follow
    And user points out treatment arms section
    # caption: The plan saved: the study will be followed against it
    When user clicks on "Save plan" button
    # caption: The plan saved: the study will be followed against it
    And user enters "MC38 guide plan" into "Name" input
    # caption: The plan saved: the study will be followed against it
    And user clicks on "Save" button in "Save plan" dialog
    # caption: The study under way: its plan and its measurements so far
    When user clicks on overview step and the stage opens
    # caption: The study under way: its plan and its measurements so far
    And user clicks on open a study strip
    # caption: The study under way: its plan and its measurements so far
    And user clicks on saved plans link
    # caption: The study under way: its plan and its measurements so far
    And user clicks on guide plan menu item
    # caption: The study under way: its plan and its measurements so far
    And user uploads "../files/example/mc38-ld4172-anti-pd1.csv" through second file zone
    # caption: Monitor: the running study against what its plan predicted
    And user clicks on still running choice
    # caption: Monitor: the running study against what its plan predicted
    And user clicks on "Open in Monitor" button in "Open a study" dialog and the stage opens
    # caption: Monitor: the running study against what its plan predicted
    Then "Needs attention" section should be visible
    # caption: Each animal against the range expected next; any of them flagged for review, with its reason
    When user clicks on "Measurements" tab
    # caption: Each animal against the range expected next; any of them flagged for review, with its reason
    And user picks "Flag this animal…" from the context menu of first sheet row
    # caption: Each animal against the range expected next; any of them flagged for review, with its reason
    And user enters "Ulcerated tumor: check at the next session" into "Reason" input
    # caption: Each animal against the range expected next; any of them flagged for review, with its reason
    And user clicks on "OK" button and PMX Studio settles
    # caption: Welfare against the humane endpoints, and every decision logged with who, when and why
    When user clicks on "Welfare" tab and PMX Studio settles
    # caption: Welfare against the humane endpoints, and every decision logged with who, when and why
    When user clicks on "Decisions" tab and PMX Studio settles
    # caption: Analyze: the finished study, every arm against the control
    When user clicks on analyze header step and the stage opens
    # caption: Analyze: the finished study, every arm against the control
    And user points out experimental groups section
    # caption: Time to the endpoint, projected past the last measurement
    And user scrolls the page to km chart
    # caption: Time to the endpoint, projected past the last measurement
    And user drags the slider of km slider to 24
    # caption: Each animal: its response and its tumor over time
    And user scrolls the page to each animal tab
    # caption: Each animal: its response and its tumor over time
    When user clicks on each animal tab and PMX Studio settles
    # caption: Each animal: its response and its tumor over time
    And user points out waterfall chart
    # caption: Each animal: its response and its tumor over time
    And user scrolls the page to twin section
    # caption: Its digital twin: the same tumor untreated; the gap is the effect of the treatment
    And user clicks on third animal row and PMX Studio settles
    # caption: Its digital twin: the same tumor untreated; the gap is the effect of the treatment
    And user points out twin effect line
    # caption: Tolerability, and the study against its plan
    And user scrolls the page to tolerability tab
    # caption: Tolerability, and the study against its plan
    When user clicks on "Tolerability" tab and PMX Studio settles
    # caption: Tolerability, and the study against its plan
    When user clicks on "Against the plan" tab and PMX Studio settles
    # caption: The whole portfolio, from design to analysis
    When user clicks on overview step and the stage opens
