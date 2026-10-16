@journey
Feature: A Hit Design campaign keeps its columns' tags and formulas when it is saved and reopened
  A campaign saves its table as CSV and, in campaign.json, each column's semantic type and type. A tag set on
  a column (".choices" makes the grid edit the column with a dropdown) and a calculated column's formula live
  in the column's tags, and the save kept neither: the reopened campaign had no dropdown, and its calculated
  column came back as plain values that no longer followed their source. The campaign is reopened from the
  app's campaign list, which reads both files from the server again.

  The campaign is made for the feature and removed at its end, with its rows in the hitdesign database
  (bindings/hit-design.ts). Its Notes column is empty, as a column the user has just added is.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Columns" with a table of:
      | Molecule | Stage  | V-iD | Score | Notes |
      | CCO      | Design |      | 10    |       |
      | c1ccccc1 | Design |      | 4     |       |
    And user opens the "Hit Design" app
    When user clicks on "BDD Columns" link
    Then the "BDD Columns" view should be current
    And the table should have 2 rows

  Scenario: A tag set on an empty column is there when the campaign is reopened
    When user picks "Column Properties..." from the context menu of the "header Notes" area of grid
    And user adds the ".choices" tag with value '["low","high"]' to "Tags" input in "Notes" dialog
    And user clicks on OK button in "Notes" dialog
    And user double-clicks on the "cell 1 of Notes" area of grid
    Then cell editor should offer the choice "low"
    When user presses Escape
    And user clicks on "Submit..." button
    And user clicks on Save button in "Submit" dialog
    Then an info balloon containing "Campaign saved successfully" should have been shown
    When user clicks on "Hit Design" link
    And user clicks on "BDD Columns" link
    Then the "BDD Columns" view should be current
    And the table should have the columns "Molecule, Stage, V-iD, Score, Notes"
    And "Notes" column should have tag ".choices" equal to '["low","high"]'
    When user double-clicks on the "cell 1 of Notes" area of grid
    Then cell editor should offer the choice "low"
    And cell editor should offer the choice "high"
    When user presses Escape
    Then no error or warning balloon should have been shown

  Scenario: A calculated column keeps its formula and follows its source when the campaign is reopened
    When user adds a calculated column "Score x2" with formula "${Score} * 2"
    Then the value of "Score x2" column in row 1 should be "20"
    When user clicks on "Submit..." button
    And user clicks on Save button in "Submit" dialog
    Then an info balloon containing "Campaign saved successfully" should have been shown
    When user clicks on "Hit Design" link
    And user clicks on "BDD Columns" link
    Then the "BDD Columns" view should be current
    And the table should have the columns "Molecule, Stage, V-iD, Score, Notes, Score x2"
    And "Score x2" column should have tag "formula" equal to "${Score} * 2"
    When user double-clicks on the "cell 1 of Score" area of grid
    And user presses Control+A in cell editor
    And user types "7" into cell editor
    And user presses Enter
    Then the value of "Score" column in row 1 should be "7"
    And the value of "Score x2" column in row 1 should be "14"
    And no error or warning balloon should have been shown
