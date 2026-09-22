@journey @realizes:chem.int.project-save-reopen-chem-state
Feature: A saved project brings back the Chem state it was saved with
  On spgi-100 a substructure search for pyridine keeps 17 of the 100 rows on a substructure card for
  Structure, and a generated Scaffold Tree with its first node checked narrows the table further.
  Saved as a project and reopened after everything is closed, the view comes back with the Scaffold
  Tree, the substructure card holding pyridine, and the same rows passing; unchecking the restored
  node widens the table back to the 17 rows the substructure search alone keeps.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens spgi dataset

  Scenario: The search and the tree narrow the table before it is saved
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then "Substructure search" dialog should be visible
    When user clicks on OK button in "Substructure search" dialog
    Then sketcher dialog should be visible
    When user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    And user hovers over Scaffold Tree viewer
    And user clicks on "Generate" icon inside Scaffold Tree viewer
    Then Scaffold Tree viewer should have finished building its tree
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 1
    And fewer than 17 rows should pass the filter
    And no errors should have been logged

  Scenario: The project comes back with the tree, the card and the rows
    When user saves the current view as project "chem-state-roundtrip"
    And user closes all views
    And user opens the "chem-state-roundtrip" project
    Then Scaffold Tree viewer should be visible
    And the "nodes" reading of Scaffold Tree viewer should be at least 1
    And the "checked nodes" reading of Scaffold Tree viewer should be 1
    And "Structure" filter card should be visible
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And the table should have 100 rows
    And fewer than 17 rows should pass the filter
    And no errors should have been logged

  Scenario: Unchecking the restored node leaves the substructure search alone
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 0
    And 17 rows should pass the filter
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccncc1"
    And no errors should have been logged
