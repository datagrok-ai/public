@journey
Feature: MPO Score counts the measured profile parameters of a row with Coverage Score
  Chem | Calculate | MPO Score... has a Coverage Score checkbox under Pareto front, off by default
  (public#4403). Checked, OK adds an int column "MPO <profile> coverage" beside the score: the number
  of the profile's properties that have a value in the row. An empty cell is not counted, whether the
  column holds plain numbers or qualified ones; a qualified value such as "<-5" is. The count does not
  depend on the score: a row the score leaves empty (Exclude row, the default for a missing value)
  still gets its count. Unchecked, OK adds no coverage column. A second checked run replaces the score
  and the coverage columns instead of adding copies.

  The table is the shipped ADME Test profile's three properties with some cells empty: row 1 has all
  three, row 2 misses Caco2, row 3 misses Caco2 and Solubility, row 4 has none, row 5 has a qualified
  Caco2, so its counts are 3, 2, 1, 0, 3. What an unchecked run does to a coverage column left by an
  earlier checked one, and a coverage option in the MPO context panel, are questions in MISSING.md;
  log10 curves and the cached calculator are asked for as package tests there.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the shipped MPO profile "ADME Test" is on the server
    And user opens a table "adme coverage" with:
      | smiles   | Caco2 | Lipophilicity | Solubility |
      | c1ccccc1 | -6    | 2             | 0          |
      | CCO      |       | 2             | 0          |
      | c1ccncc1 |       | 3             |            |
      | CC(=O)O  |       |               |            |
      | CCN      | <-5   | 2             | 0          |

  Scenario: Coverage Score is off by default and a run without it adds no coverage column
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    Then "MPO Score" dialog should be visible
    And "Coverage Score" input in "MPO Score" dialog should be below "Pareto front" input in "MPO Score" dialog
    And "Coverage Score" input in "MPO Score" dialog should have value "false"
    When user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And a new column "MPO ADME Test" should have been added
    And 1 new column should have been added
    And no errors should have been logged

  Scenario: Checked, the run adds the coverage column with the count of measured parameters per row
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user checks "Coverage Score" input in "MPO Score" dialog
    And user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And no error or warning balloon should have been shown
    And a new column "MPO ADME Test coverage" should have been added
    And "MPO ADME Test coverage" column should have type "int"
    And the "MPO ADME Test coverage" cell of row 1 should be displayed as "3"
    And the "MPO ADME Test coverage" cell of row 2 should be displayed as "2"
    And the "MPO ADME Test coverage" cell of row 3 should be displayed as "1"
    And the "MPO ADME Test coverage" cell of row 4 should be displayed as "0"
    And the "MPO ADME Test coverage" cell of row 5 should be displayed as "3"
    And the "MPO ADME Test" cell of row 2 should be displayed as ""
    And no errors should have been logged

  Scenario: A second run with Coverage Score replaces both columns instead of adding copies
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user checks "Coverage Score" input in "MPO Score" dialog
    And user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And the table should have 6 columns
    And the "MPO ADME Test coverage" cell of row 4 should be displayed as "0"
    And no errors should have been logged
