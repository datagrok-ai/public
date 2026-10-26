@realizes:chem.cp.scaffold-tree-add-filter
Feature: Scaffold Tree regressions from closed issues
  Closed issues of the Scaffold Tree viewer and its filter card on spgi-100 (100 rows, molecule column
  Structure), each a scenario of its own on a fresh table. Every tree is sketched by hand ("Sketch
  scaffolds manually"), so no scenario runs the Python tree generation. Closing a viewer whose scaffold
  is colored drops the "Structure colors" column and the highlight tag on Structure (public#3408). The
  filter panel's reset unchecks the checked scaffold (public#2459, public#2154). Select | Selection to
  Filter holds its filter while a scaffold tree is in the view (public#2964). Order or Hide Columns lists
  the same columns after a tree with a checked scaffold is added (public#2315). A child added to a
  scaffold drawn with an explicit H opens with no "not a superstructure" warning, while a molecule that
  is not one shows it (public#2626). Editing one scaffold keeps another one checked (public#2135, here on
  the viewer: the filter card's readings are not reachable from Chem's bindings, see MISSING.md). A
  property changed from the context panel keeps the scatter plot's properties there while a scaffold
  tree filter holding a scaffold is on the panel (public#2322).

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens spgi dataset

  Scenario: Closing a viewer with a colored scaffold drops the colors column and the highlight (public#3408)
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    Then Scaffold Tree viewer should be visible
    When user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 1
    When user hovers over the "node 1" area of Scaffold Tree viewer
    And user clicks on the "color icon of node 1" area of Scaffold Tree viewer
    Then the "colored nodes" reading of Scaffold Tree viewer should be 1
    And the table should have a column "Structure colors"
    When user clicks on close icon of Scaffold Tree viewer
    Then the open tableview should have 0 Scaffold Tree viewers
    And the table should not have a column "Structure colors"
    And "Structure" column should have tag ".chem-scaffold-tree-highlight" equal to "[]"
    And no errors should have been logged

  Scenario: The filter panel's reset unchecks the checked scaffold (public#2459, public#2154)
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    And user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 1
    And the "hits of node 1" reading of Scaffold Tree viewer should be 17
    When user clicks on filter icon in toolbar
    And user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 1
    And 17 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the "checked nodes" reading of Scaffold Tree viewer should be 0
    And the "checked of node 1" reading of Scaffold Tree viewer should be "false"
    And no errors should have been logged

  Scenario: Selection to Filter holds its filter while a scaffold tree is in the view (public#2964)
    When user clicks on filter icon in toolbar
    And user picks "Chem > Analyze > Scaffold Tree" from the top menu
    And user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccccc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 1
    When user selects the first 1 rows
    And user picks "Select > Selection to Filter" from the top menu
    Then 1 row should pass the filter
    And the "checked nodes" reading of Scaffold Tree viewer should be 0
    And the "nodes" reading of Scaffold Tree viewer should be 1
    And 1 row should pass the filter
    And no errors should have been logged

  Scenario: Order or Hide Columns lists the same columns after a tree with a checked scaffold is added (public#2315)
    When user picks "Order or Hide Columns..." from the context menu of the "header Id" area of grid
    Then Order or Hide Columns dialog should be visible
    And the "rows shown" reading of Grid viewer in Order or Hide Columns dialog should be 88
    When user remembers the "rows shown" reading of Grid viewer in Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    And user picks "Chem > Analyze > Scaffold Tree" from the top menu
    And user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "hits of node 1" reading of Scaffold Tree viewer should be 17
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then 17 rows should pass the filter
    When user picks "Order or Hide Columns..." from the context menu of the "header Id" area of grid
    Then Order or Hide Columns dialog should be visible
    And the "rows shown" reading of Grid viewer in Order or Hide Columns dialog should be as remembered
    When user clicks on CLOSE button in Order or Hide Columns dialog
    Then no errors should have been logged

  Scenario: A child of a scaffold drawn with an explicit H opens with no superstructure warning (public#2626)
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    And user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "[H]c1ccccc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 1
    When user hovers over the "node 1" area of Scaffold Tree viewer
    And user clicks on the "add icon of node 1" area of Scaffold Tree viewer
    Then sketcher dialog should be visible
    And the sketcher in sketcher dialog should hold the molecule "[H]c1ccccc1"
    And sketcher dialog should not contain text "not a superstructure"
    When user clears molecule input of sketcher dialog
    And user types "C1CCNCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then sketcher dialog should contain text "The edited molecule is not a superstructure of its parent"
    When user clears molecule input of sketcher dialog
    And user types "Cc1ccccc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then sketcher dialog should not contain text "not a superstructure"
    When user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 2
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Editing one scaffold keeps another one checked (public#2135)
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    And user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 1
    When user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccccc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 2
    And the "scaffold of node 2" reading of Scaffold Tree viewer should be the molecule "c1ccccc1"
    And the "hits of node 2" reading of Scaffold Tree viewer should be 32
    When user clicks on the "checkbox of node 2" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 1
    And 32 rows should pass the filter
    And user hovers over the "node 1" area of Scaffold Tree viewer
    And user clicks on the "edit icon of node 1" area of Scaffold Tree viewer
    Then sketcher dialog should be visible
    When user clears molecule input of sketcher dialog
    And user types "c1ccc2ncccc2c1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Save" button in sketcher dialog
    Then the "scaffold of node 1" reading of Scaffold Tree viewer should be the molecule "c1ccc2ncccc2c1"
    And the "checked of node 2" reading of Scaffold Tree viewer should be "true"
    And the "checked nodes" reading of Scaffold Tree viewer should be 1
    And 32 rows should pass the filter
    And no errors should have been logged

  Scenario: A property changed with a scaffold tree filter on the panel keeps the viewer's properties in the context panel (public#2322)
    When user clicks on filter icon in toolbar
    And user picks "Add Filter | Scaffold Tree Filter..." from the viewer menu of filter panel
    Then "Select columns..." dialog should be visible
    And the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "Structure"
    When user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user hovers over filter panel
    And user clicks on "Sketch scaffolds manually" icon inside filter panel
    And user types "c1ccccc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then sketcher dialog should be hidden
    Given user adds a scatter plot viewer
    When user clicks on settings icon of scatter plot viewer
    And user checks "Show Filtered Out Points" property in context panel
    Then "showFilteredOutPoints" property of scatter plot viewer should be "true"
    And "Show Filtered Out Points" property in context panel should be visible
    When user unchecks "Show Filtered Out Points" property in context panel
    Then "showFilteredOutPoints" property of scatter plot viewer should be "false"
    And no errors should have been logged
