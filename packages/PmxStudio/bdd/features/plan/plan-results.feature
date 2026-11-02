@journey
Feature: Plan's results
  The figures' tabs one by one: the volume axis, the calendar's weekdays and randomization day, the trade-offs, the
  protocol ready to paste, saving and downloading the plan, and moving between the stages.

  Background:
    Given user is logged in
    And no plan named "BDD plan" is saved
    And user opens PMX Studio
    Then the study should have finished running

  Scenario: Design tab
    Then "Experimental groups" section should contain text "Treatment 1"
    And "Predicted tumor growth" section should be visible
    And visible volume axis should have the value "Log"
    When user clicks on linear axis choice
    Then visible volume axis should have the value "Linear"
    When user clicks on log axis choice
    Then visible volume axis should have the value "Log"

  Scenario: Calendar tab, measurement days
    When user clicks on "Calendar" tab
    Then "Study calendar" section should be visible
    And calendar summary should contain text "9 measurement days"
    When user clicks on tuesday header
    Then the study should have finished running
    And calendar summary should contain text "12 measurement days"
    When user clicks on tuesday header
    Then calendar summary should contain text "9 measurement days"

  Scenario: Calendar tab, randomization day
    When user selects "Wednesday" in randomization day picker
    Then calendar summary should contain text "Randomize on a"
    And randomization day picker should have the value "Wednesday"
    When user selects "Monday" in randomization day picker
    Then randomization day picker should have the value "Monday"
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "days 0, 2, 4, 7"
    And randomization day picker should have the value "Monday"

  Scenario: Calendar tab, two changes without waiting are both planned
    When user clicks on "Calendar" tab
    And user clicks on tuesday header
    And user clicks on thursday header
    And user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "Mon, Tue, Wed, Thu, Fri"
    When user clicks on "Calendar" tab
    And user clicks on tuesday header
    And user clicks on thursday header
    And user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "on Mon, Wed, Fri"

  Scenario: A form edited during a replan stays out of date
    When user clicks on "Calendar" tab
    And user selects "Tuesday" in randomization day picker
    And user enters "1400" into "Endpoint volume" input
    Then the study should have finished running
    And randomization day picker should have the value "Tuesday"
    And stale note should be visible
    When user enters "1500" into "Endpoint volume" input
    And user selects "Monday" in randomization day picker
    And user clicks on form run button
    Then the study should have finished running
    And stale note should be hidden

  Scenario: Trade-offs tab
    When user clicks on "Trade-offs" tab
    Then "Power by group size" section should be visible
    And "Power by end of study" section should contain text "Ending after day"
    And "Measurement schedule" section should contain text "this plan"
    When user clicks on "Design" tab
    And user enters "1400" into "Endpoint volume" input
    And user clicks on form run button
    Then the study should have finished running
    When user enters "1500" into "Endpoint volume" input
    And user clicks on "Trade-offs" tab
    Then "Measurement schedule" section should contain text "Run to compare"
    When user clicks on form run button
    Then the study should have finished running
    And "Measurement schedule" section should contain text "this plan"

  Scenario: Protocol tab
    When user clicks on "Protocol" tab
    Then "Study protocol" section should contain text "MC38"
    And "Study protocol" section should contain text "to implant"
    When user clicks on "Copy" button
    Then an info balloon containing "Study protocol copied" should have been shown
    And the clipboard should contain text "Study protocol (Plan MC38"

  Scenario: Saving a plan, and replacing it on a second press
    When user clicks on "Save plan" button
    Then "Save plan" dialog should be visible
    When user enters "BDD plan" into "Name" input in "Save plan" dialog
    And user clicks on "Save" button in "Save plan" dialog
    Then an info balloon containing "Plan saved" should have been shown
    And a plan named "BDD plan" should be saved
    And "Save plan" dialog should be absent
    When user clicks on "Save plan" button
    And user enters "BDD plan" into "Name" input in "Save plan" dialog
    And user clicks on "Save" button in "Save plan" dialog
    Then "Save plan" dialog should contain text "already saved"
    When user clicks on "Save" button in "Save plan" dialog
    Then "Save plan" dialog should be absent

  Scenario: Downloading the plan
    Given user watches downloads
    When user downloads a file through "Download" button
    Then a file matching "\.pmxplan$" should have been downloaded
    And the downloaded file should contain "pmx-plan"

  Scenario: A plan out of date is not saved
    When user enters "1200" into "Endpoint volume" input
    And user clicks on "Save plan" button
    Then a warning balloon containing "run to update the plan first" should have been shown
    And "Save plan" dialog should be absent
    When user enters "1500" into "Endpoint volume" input
    And user clicks on form run button
    Then the study should have finished running

  Scenario: Moving between the stages
    When user clicks on "3 Analyze" button
    Then analyze start should be visible
    When user clicks on "1 Plan" button
    Then "Study protocol" section should be visible
