@journey @realizes:chem.cp.scaffold-tree-add-filter
Feature: Scaffold Tree colors, blocked generation and two tables
  A Scaffold Tree generated on spgi-100: coloring its first node adds
  the hidden "Structure colors" column and names it in the viewer's reading, and a scatter plot
  colored by that column takes the scaffold's colors. The magic wand is blocked, with its reason, on
  a table with no molecule column and on a molecule column of 500 categories or more. With two tables
  open the menu binds a new viewer to the active one. Removing another scaffold keeps the colors
  column and the plot's coloring; removing the last colored one drops the column, and the plot
  goes back to no coloring.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens spgi dataset
    And user picks "Chem > Analyze > Scaffold Tree" from the top menu

  Scenario: Coloring a scaffold adds the colors column
    When user hovers over Scaffold Tree viewer
    And user clicks on "Generate" icon inside Scaffold Tree viewer
    Then Scaffold Tree viewer should have finished building its tree
    And the table should not have a column "Structure colors"
    When user hovers over the "node 1" area of Scaffold Tree viewer
    And user clicks on the "color icon of node 1" area of Scaffold Tree viewer
    Then the "colored nodes" reading of Scaffold Tree viewer should be 1
    And the table should have a column "Structure colors"
    And the "colors column" reading of Scaffold Tree viewer should be "Structure colors"
    And "Structure colors" column should have at least 1 distinct values
    And "Structure colors" column should be color-coded categorically
    And no errors should have been logged

  Scenario: A scatter plot colored by the colors column takes the scaffold's colors
    Given user adds a scatter plot viewer
    When user picks "Structure colors" in the "color" column selector of scatter plot viewer
    Then "colorColumnName" property of scatter plot viewer should be "Structure colors"
    And the legend of scatter plot viewer should list 2 items
    And no errors should have been logged

  Scenario: The magic wand is blocked on a table with no molecule column
    Given user opens demog dataset
    And user adds Scaffold Tree viewer
    Then Scaffold Tree viewer should be visible
    And the "generate blocked reason" reading of Scaffold Tree viewer should be "There is no molecule column in the table"
    And the "message" reading of Scaffold Tree viewer should include the text "No molecule column found"
    When user hovers over Scaffold Tree viewer
    Then "Generate" icon inside Scaffold Tree viewer should be disabled
    And the "nodes" reading of Scaffold Tree viewer should be 0
    And no errors should have been logged

  Scenario: The magic wand is blocked past 500 structure categories
    Given user opens mol1K dataset
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    Then Scaffold Tree viewer should be visible
    And the "generate blocked reason" reading of Scaffold Tree viewer should be "The number of molecules exceeds the limit of 500"
    When user hovers over Scaffold Tree viewer
    Then "Generate" icon inside Scaffold Tree viewer should be disabled
    And the "nodes" reading of Scaffold Tree viewer should be 0
    And no errors should have been logged

  Scenario: With two tables open the menu binds the viewer to the active one (github-3004)
    Given user opens spgi dataset keeping the first 20 rows as "tableA"
    And user opens spgi dataset keeping the first 30 rows as "tableB"
    Then the "tableB" view should be current
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    Then Scaffold Tree viewer should be bound to table "tableB"
    And the open tableview should have 1 Scaffold Tree viewer
    When user switches to the "tableA" table view
    Then the open tableview should have 0 Scaffold Tree viewers
    And no errors should have been logged

  Scenario: Removing an uncolored scaffold keeps the colors column and the plot colored by it
    Given user switches to the "spgi-100" table view
    Then the "colored nodes" reading of Scaffold Tree viewer should be 1
    When user remembers the "nodes" reading of Scaffold Tree viewer
    And user hovers over the "node 2" area of Scaffold Tree viewer
    And user clicks on the "remove icon of node 2" area of Scaffold Tree viewer
    And user clicks on "Yes" button in "Remove scaffold" dialog
    Then the "nodes" reading of Scaffold Tree viewer should not be as remembered
    And the "colored nodes" reading of Scaffold Tree viewer should be 1
    And the table should have a column "Structure colors"
    And "colorColumnName" property of scatter plot viewer should be "Structure colors"
    And no errors should have been logged

  Scenario: Removing the last colored scaffold drops the colors column and the plot's coloring
    Then the "color of node 1" reading of Scaffold Tree viewer should not be ""
    When user hovers over the "node 1" area of Scaffold Tree viewer
    And user clicks on the "remove icon of node 1" area of Scaffold Tree viewer
    And user clicks on "Yes" button in "Remove scaffold" dialog
    Then the "colored nodes" reading of Scaffold Tree viewer should be 0
    And the table should not have a column "Structure colors"
    And "colorColumnName" property of scatter plot viewer should be ""
    And scatter plot viewer should be painted
    And no errors should have been logged
