@journey @realizes:filters.cp.chem-and-bio-filters
Feature: A structure drawn in Ketcher reaches the substructure card as drawn
  Chem | Search | Substructure Search... on smiles (1000 rows, molecule column canonical_smiles) adds a
  substructure card and opens its sketcher, here Ketcher. The benzene template placed on the canvas
  keeps the 924 rows with a benzene ring. Two ways the drawing got lost: with Filter as you draw
  cleared, OK right after the stroke read the molecule before Ketcher's asynchronous export had it,
  and the dialog's closing dropped that export, so the card stayed empty; and while the pointer is
  over the canvas, Ketcher holds the template's floating preview in its structure, so an export taken
  then read two rings (666 rows). The pointer drifts on over the canvas after every stroke here.

  Background:
    Given user is logged in
    And the molecule sketcher is "Ketcher"
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: With Filter as you draw cleared, OK right after the stroke filters by it
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then sketcher dialog should be visible
    When user unchecks "Filter as you draw" input in sketcher dialog
    And user places the benzene template on the Ketcher canvas
    And user clicks on OK button in sketcher dialog
    Then 924 rows should pass the filter
    And the "structure of canonical_smiles" reading of filter panel should be "c1ccccc1"
    And "Sketch" text in "canonical_smiles" filter card should not be visible
    And no errors should have been logged

  Scenario: With Filter as you draw checked, the rows follow the stroke while the pointer rests on the canvas
    When user hovers over "canonical_smiles" filter card
    And user clicks on close of "canonical_smiles" filter card
    Then all rows should pass the filter
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then sketcher dialog should be visible
    When user checks "Filter as you draw" input in sketcher dialog
    And user places the benzene template on the Ketcher canvas
    Then 924 rows should pass the filter
    And the "structure of canonical_smiles" reading of filter panel should be "c1ccccc1"
    When user clicks on OK button in sketcher dialog
    Then 924 rows should pass the filter
    And the "structure of canonical_smiles" reading of filter panel should be "c1ccccc1"
    And no errors should have been logged
