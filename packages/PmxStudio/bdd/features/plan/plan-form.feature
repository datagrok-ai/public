@journey
Feature: Plan's declaration form
  Every field of the form, one by one: it takes a value, the plan is marked out of date, Run brings the change into the
  figures, and a value out of range holds Run back. The numbers themselves are the package tests'.

  Background:
    Given user is logged in
    And user opens PMX Studio
    Then the study should have finished running
    And "Experimental groups" section should contain text "Control"

  Scenario: A loaded plan is not out of date and Run is on
    Then stale note should be hidden
    And form run button should be enabled
    And "Cell line" input should have the value "MC38"
    And "Endpoint volume" input should have the value "1500"
    And "Target power" input should have the value "80"

  Scenario: Cell line
    When user selects "CT26" in "Cell line" input
    Then stale note should be visible
    When user clicks on form run button
    Then the study should have finished running
    And stale note should be hidden
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "CT26"
    When user clicks on "Design" tab
    And user selects "B16F10" in "Cell line" input
    And user clicks on form run button
    Then the study should have finished running
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "B16F10"
    When user clicks on "Design" tab
    And user selects "MC38" in "Cell line" input
    And user clicks on form run button
    Then the study should have finished running

  Scenario: Endpoint volume
    When user remembers the "End of study" figure
    And user enters "800" into "Endpoint volume" input
    Then stale note should be visible
    When user clicks on form run button
    Then the study should have finished running
    And the "End of study" figure should be lower than remembered
    When user enters "100" into "Endpoint volume" input
    Then form run button should be disabled
    When user enters "5000" into "Endpoint volume" input
    Then form run button should be disabled
    When user clears "Endpoint volume" input
    Then form run button should be disabled
    When user enters "1500" into "Endpoint volume" input
    Then form run button should be enabled
    When user clicks on form run button
    Then the study should have finished running

  Scenario: Target power
    When user remembers the "Analyzable" figure
    And user enters "90" into "Target power" input
    And user clicks on form run button
    Then the study should have finished running
    And the "Analyzable" figure should be higher than remembered
    When user enters "99" into "Target power" input
    And user clicks on form run button
    Then the study should have finished running
    And study error should be hidden
    And the "Analyzable" figure should read "—"
    And "Experimental groups" section should contain text "do not reach 99% power"
    When user clicks on "Trade-offs" tab
    Then "Power by end of study" section should contain text "No group size"
    And "Measurement schedule" section should contain text "No group size"
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "Not sized"
    When user clicks on "Design" tab
    When user enters "10" into "Target power" input
    Then form run button should be disabled
    When user enters "80" into "Target power" input
    And user clicks on form run button
    Then the study should have finished running

  Scenario: Animals per arm
    When user enters "12" into "Animals per arm" input
    And user clicks on form run button
    Then the study should have finished running
    And the "Analyzable" figure should read "24"
    And "Experimental groups" section should contain text "Smallest detectable effect"
    And the smallest detectable effects should have arrived
    When user enters "1" into "Animals per arm" input
    Then form run button should be disabled
    When user clears "Animals per arm" input
    Then form run button should be enabled
    When user clicks on form run button
    Then the study should have finished running
    And the "Analyzable" figure should read "16"

  Scenario: End of study
    When user enters "25" into "End of study" input
    And user clicks on form run button
    Then the study should have finished running
    And the "End of study" figure should read "Day 25"
    When user enters "3" into "End of study" input
    Then form run button should be disabled
    When user clears "End of study" input
    And user clicks on form run button
    Then the study should have finished running
    And the "End of study" figure should read "Day 18"

  Scenario: Advanced fields
    Then "Significance (α)" input should be hidden
    When user clicks on advanced header
    Then "Significance (α)" input should be visible
    And "Lost before randomization" input should be visible
    And "Lost after randomization" input should be visible
    And "Extra animals for losses" input should be visible
    And "Confidence" input should be visible
    And "Volume at randomization" input should be visible
    And "Humane limit" input should be visible
    And "Measurement days" input should be hidden
    And "Randomization day" input should be hidden
    And "Compare measurement schedules" input should be hidden

  Scenario: Losses
    When user remembers the "Animals to implant" figure
    And user remembers the "To randomize" figure
    And user enters "50" into "Lost before randomization" input
    And user clicks on form run button
    Then the study should have finished running
    And the "Animals to implant" figure should be higher than remembered
    And the "To randomize" figure should be as remembered
    When user enters "20" into "Lost before randomization" input
    And user enters "30" into "Lost after randomization" input
    And user clicks on form run button
    Then the study should have finished running
    And the "To randomize" figure should be higher than remembered
    When user enters "10" into "Lost after randomization" input
    And user selects "Average loss only" in "Extra animals for losses" input
    And user clicks on form run button
    Then the study should have finished running
    And the "Animals to implant" figure should be lower than remembered
    When user selects "Enough for every group" in "Extra animals for losses" input
    And user enters "99" into "Confidence" input
    And user clicks on form run button
    Then the study should have finished running
    And the "Animals to implant" figure should be higher than remembered
    When user enters "90" into "Confidence" input
    And user enters "120" into "Lost before randomization" input
    Then form run button should be disabled
    When user enters "20" into "Lost before randomization" input
    And user clicks on form run button
    Then the study should have finished running

  Scenario: Volume at randomization, significance and humane limit
    When user remembers the "End of study" figure
    And user enters "300" into "Volume at randomization" input
    And user clicks on form run button
    Then the study should have finished running
    And the "End of study" figure should be lower than remembered
    When user enters "100" into "Volume at randomization" input
    And user enters "0.01" into "Significance (α)" input
    And user clicks on form run button
    Then the study should have finished running
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "α = 0.01"
    When user clicks on "Design" tab
    And user enters "0.05" into "Significance (α)" input
    And user enters "1000" into "Humane limit" input
    And user clicks on form run button
    Then the study should have finished running
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "1000 mm³"
    When user clicks on "Design" tab
    And user enters "2000" into "Humane limit" input
    And user clicks on form run button
    Then the study should have finished running
