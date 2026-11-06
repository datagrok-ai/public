@guide
Feature: PMX Studio in two minutes
  A teaser: the answer to "what does PMX Studio do for a preclinical efficacy study?". A lab's portfolio on the
  Overview; one study, MC38, at its three moments: its plan, followed halfway through (an animal flagged, the forecasts
  recomputed) and analyzed with each animal against its digital twin and followed past its last measurement.

  Scenario: One study from its plan to its analysis
    Given the guide may take up to 8 minutes
    And user is logged in
    And simple mode is off
    And the shell is in presentation mode without tooltips
    # caption: PMX Studio: the lab's studies, planned, running and finished
    And user opens PMX Studio
    Then the studies should be listed
    # caption: PMX Studio: the lab's studies, planned, running and finished
    And user waits 3 seconds
    # caption: The plan of an MC38 study, sized from the atlas of historical controls
    When user clicks on demo plan row
    Then the plan stage should be open
    # caption: The plan of an MC38 study, sized from the atlas of historical controls
    And user waits 2 seconds
    # caption: Change an assumption: the design follows at once
    When user enters "1200" into "Endpoint volume" input
    Then the study should have finished running
    # caption: Change an assumption: the design follows at once
    And user waits 2 seconds
    # caption: The same study, halfway through
    When user clicks on overview step
    # caption: The same study, halfway through
    And user clicks on demo followed row
    Then the monitor stage should be open
    # caption: Monitor: who needs attention today, and when the study should end
    And user waits 3 seconds
    # caption: Flag an animal from its row
    When user picks "Flag this animal…" from the context menu of first attention row
    # caption: With its reason: it goes in the study's log
    And user enters "Ulcerated tumor: check at the next session" into "Reason" input
    # caption: With its reason: it goes in the study's log
    And user clicks on "OK" button
    Then the study should have finished running
    # caption: The forecasts are recomputed with it
    And user waits 2 seconds
    # caption: The study with all its data
    When user clicks on overview step
    # caption: The study with all its data
    And user clicks on demo finished row
    Then the analyze stage should be open
    # caption: Analyze: every arm against the control
    And user waits 2 seconds
    # caption: Each animal against its digital twin, untreated
    When user clicks on each animal tab
    Then the study should have finished running
    # caption: Each animal against its digital twin, untreated
    And user waits 3 seconds
    # caption: Survival to the endpoint
    When user clicks on results tab
    # caption: Survival to the endpoint
    And user scrolls km chart into view
    # caption: Follow the study past its last measurement
    And user drags the slider of km slider to 40
    # caption: A projection, marked as such
    Then "Time to endpoint" section should contain text "projected"
    # caption: A projection, marked as such
    And user waits 2 seconds
    # caption: Next: combination therapies. A, B and A+B: synergy or antagonism
    When user points at "Analyze > Combination synergy…" in the top menu
    # caption: Back on the Overview: planned, followed and finished, in one place
    When user clicks on overview step
    Then the studies should be listed
    # caption: Back on the Overview: planned, followed and finished, in one place
    And user waits 3 seconds
