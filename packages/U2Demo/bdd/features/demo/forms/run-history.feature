@demo @realizes:u2.function-input @realizes:u2.func-call-history-browser
Feature: Run history
  A function input builds a funcForm with the standard run button and history icon; a run is
  saved on the server and comes back in the history popup.

  Scenario: Running a function and opening its history
    Given user opens the "Run history" demo page
    When user selects "Abs" in Function input
    Then function form should be visible
    And run button of function form should be visible
    When user enters "3" into x input in function form
    And user clicks on run button of function form
    And user clicks on history icon of function form
    Then history browser should be visible
