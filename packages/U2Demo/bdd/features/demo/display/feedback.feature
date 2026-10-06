@demo @realizes:u2.progress @realizes:u2.notify @realizes:u2.tour
Feature: Progress, notifications and the guided tour

  Background:
    Given user opens the "Feedback" demo page

  Scenario: Progress runs to completion and resets
    Then "Indexing compounds" progress bar should contain text "0%"
    When user clicks on Run button
    Then "Indexing compounds" progress bar should contain text "100%"
    And value of progress readout should have text "1.00"
    When user clicks on Reset button
    Then value of progress readout should have text "0.00"

  Scenario: Notifications stack and close
    When user clicks on Info button
    Then notification should contain text "Indexed 1,204 compounds."
    When user clicks on Warning button
    And user clicks on Error button
    Then "Two structures could not be parsed." notification should be visible
    And "Connection refused." notification should be visible
    When user closes "Connection refused." notification
    Then "Connection refused." notification should be hidden
    When user clicks on "Close all" button
    Then "Two structures could not be parsed." notification should be hidden

  Scenario: The tour runs to its end and skips a missing target
    When user clicks on "Start tour" button
    Then tour should be visible
    And "1 / 4" text should be visible
    When user clicks on NEXT button
    And user clicks on NEXT button
    And user clicks on NEXT button
    Then tour should be hidden
    And value of tour readout should have text "done"
    When user clicks on "Start tour" button
    And user presses Escape
    Then value of tour readout should have text "skipped"
