@edge-cases @realizes:views.queries
Feature: Query parameters edge cases: the run dialog refuses what it cannot run
  Guards for the dialog a parameterized query opens before it runs. Every parameter is required
  unless it is nullable, and the input is validated before the query runs: "you will get an error if
  a required parameter is not specified" (help `datagrok/concepts/functions/func-params-annotation.md`,
  "Initial values and optional parameters"). The validation also checks the min-max conditions and a
  parameter's own `validator` expression (the same page, "Validation"); a choice parameter offers its
  choices (the same page, "Choices"). The dialog shows a refused value by disabling OK, and a valid
  value enables it again. A nullable parameter may be left empty.

  The query is the feature's own, on the platform's database (System:Datagrok, its entity_types
  table, where the entity type "Project" is one row on every stand), saved through the JS API with
  the run's time in its name and deleted when the feature ends. The dialog is opened through Run in
  the Browse tree; nothing is run except in the nullable scenario. That an empty integer and an
  empty string are refused alike is core#229 (open, the owner's decision is pending): see MISSING.md.

  Background:
    Given user is logged in
    And the browse panel is open
    And a query "BDDQParamEdge{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      --input: int lim = 5 {min: 1; max: 10}
      --input: string kind = "a" {choices: ["a", "b"]}
      --input: string code = "1234" {validator: startsWith(value, "12")}
      --input: string note = "x" {nullable: true}
      select name from entity_types where name = @typeName limit @lim
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree
    Then "BDDQParamEdge{time}" dialog should be visible
    And OK button in "BDDQParamEdge{time}" dialog should be enabled

  Scenario: A cleared required parameter disables OK, and a value enables it again
    When user clears "Type Name" input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be disabled
    When user enters "Project" into "Type Name" input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clicks on CANCEL button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And no errors should have been logged

  Scenario: An integer above its maximum disables OK, and one inside the bounds enables it
    When user enters "11" into Lim input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be disabled
    When user enters "10" into Lim input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user enters "0" into Lim input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be disabled
    When user clicks on CANCEL button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And no errors should have been logged

  Scenario: Text in an integer parameter disables OK
    When user enters "abc" into Lim input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be disabled
    When user enters "3" into Lim input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clicks on CANCEL button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And no errors should have been logged

  Scenario: A value its validator rejects disables OK
    When user enters "99" into Code input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be disabled
    When user enters "1299" into Code input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clicks on CANCEL button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And no errors should have been logged

  Scenario: A choice parameter offers its choices only
    Then Kind input in "BDDQParamEdge{time}" dialog should offer the choice "b"
    And Kind input in "BDDQParamEdge{time}" dialog should not offer the choice "c"
    When user clicks on CANCEL button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And no errors should have been logged

  Scenario: A nullable parameter left empty keeps OK enabled, and the query runs
    When user clears Note input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clicks on OK button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And the current view should be a TableView view
    And the table should have 1 row
    And no errors should have been logged
    And no error or warning balloon should have been shown
