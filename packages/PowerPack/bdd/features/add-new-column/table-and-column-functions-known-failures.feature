@realizes:powerpack.dialogs.add-new-column @realizes:GROK-20931
Feature: Add New Column: what the table and column functions still get wrong
  The Add New Column dialog on a three-row table written in the feature, x = 1, 2, 3. A string argument
  naming a table the workspace has not got is refused the way ${nope} is: the dialog names it and OK is
  disabled (public#4404). On auto,
  the Type choice of Avg($[x]) names the type Avg returns, double, which is the type OK adds (public#4406).
  Found in the manual QA round of 26 Oct 2026 on GROK-20931; each scenario fails at the claim its issue
  describes until the issue is fixed.

  Not translated: that a formula naming its column as a string (Column("x"), Value("x", 0)) recalculates
  after an edit (core#1128) is the ApiTests test "DataFrame: Calculated columns: A formula naming its
  column as a string recalculates after an edit"; the 20-row preview of a table-level function
  (public#4405) is a design question, not a claim yet.

  Background:
    Given user is logged in
    And user opens a table "funcs" with:
      | x |
      | 1 |
      | 2 |
      | 3 |

  @known-failure:27 @realizes:public#4404
  Scenario: A string argument naming a table the workspace has not got is refused
    When user clicks on "Add New Column..." icon
    And user types "RowCount(\"nope\")" into formula editor
    Then formula error should contain text "nope"
    And OK button in "Add New Column" dialog should be disabled
    And no errors should have been logged

  @known-failure:35 @realizes:public#4406
  Scenario: On auto, the Type choice of Avg names the type of the column OK adds
    When user clicks on "Add New Column..." icon
    And user types "Avg($[x])" into formula editor
    Then column type input should have value "auto (double)"
    And no errors should have been logged
