@guide
Feature: PMX Studio
  A teaser: the answer to "what does PMX Studio do for a preclinical efficacy study?". The lab's portfolio on the
  Overview; an MC38 study planned, followed halfway through (an animal flagged) and analyzed, each animal against its
  digital twin and the study followed past its last measurement; moving between the stages from the header.

  Scenario: An efficacy study from its design to its analysis
    Given the guide may take up to 8 minutes
    And user is logged in
    And simple mode is off
    And the shell is in presentation mode without tooltips
    And user opens PMX Studio on its Overview
    # caption: Every tumor efficacy study of the lab, from its design to its analysis
    Then the lab's studies should be on the Overview
    # caption: Plan: group sizes, duration and power, predicted from historical control studies
    When user clicks on demo plan row and the stage opens
    # caption: Every assumption is live: lower the endpoint volume and the design follows
    And user enters "1200" into "Endpoint volume" input and PMX Studio settles
    # caption: Monitor: the running study against what its plan predicted
    When user clicks on monitor step and the stage opens
    # caption: Monitor: the running study against what its plan predicted
    And user clicks on "Follow an example study" button and PMX Studio settles
    # caption: Halfway through: each group against its forecast, and when the study should end
    Then "Needs attention" section should be visible
    # caption: For every animal, the range its next measurement should fall in
    When user clicks on "Measurements" tab
    # caption: Any animal can be flagged for review, with its reason
    And user picks "Flag this animal…" from the context menu of first sheet row
    # caption: Any animal can be flagged for review, with its reason
    And user enters "Ulcerated tumor: check at the next session" into "Reason" input
    # caption: Flags and exclusions are logged with who, when and why, and travel with the study
    And user clicks on "OK" button and PMX Studio settles
    # caption: Analyze: the finished study, every arm against the control
    When user clicks on analyze header step and the stage opens
    # caption: Analyze: the finished study, every arm against the control
    And user clicks on "Analyze an example study" button and PMX Studio settles
    # caption: Each animal against its digital twin: the same tumor, predicted untreated
    When user clicks on each animal tab and PMX Studio settles
    # caption: Time to the endpoint, as measured
    When user clicks on results tab
    # caption: Time to the endpoint, as measured
    And user scrolls km chart into view
    # caption: … and projected past the last measurement, marked as a projection
    And user drags the slider of km slider to 40
    # caption: Coming next: combination therapies, A + B against A and B, synergy or antagonism
    When user points at "Analyze > Combination synergy…" in the top menu
    # caption: The whole portfolio in one place: planned, running and finished
    When user clicks on overview step and the stage opens
