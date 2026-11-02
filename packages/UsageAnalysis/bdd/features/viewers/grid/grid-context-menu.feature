@viewers @realizes:viewers.grid
Feature: The context menu acts on the right-clicked cell
  A right click on a cell makes that cell current, as a left click does, wherever the grid is
  scrolled to and whichever row was current before. The Current Value group of the menu then acts
  on the cell under the pointer: Chem's Copy as SMILES copies that molecule, Helm's Edit Helm...
  opens that peptide. Reported 2026-09-22: with a current row set and the grid scrolled away from
  it, a right click made the clicked cell current for a moment, then the grid jumped back to the old
  current row and scrolled it into view, and the Current Value actions took whichever cell the
  pointer landed on after the jump.

  Background:
    Given user is logged in

  Scenario: A right click below the current row makes the clicked row current and keeps the scroll
    Given user opens demog-1000 dataset
    Then grid should show 1000 rows
    When user clicks on the "cell 5 of AGE" area of grid
    Then row 5 should be current
    When user scrolls the mouse wheel down 60 times over the "cell 5 of AGE" area of grid
    Then grid should have a "cell 1000 of AGE" area
    And grid should not have a "cell 5 of AGE" area
    When user right-clicks on the "cell 1000 of AGE" area of grid
    Then the open menu should list "Current Column"
    And row 1000 should be current
    And the "current row" reading of grid should be 1000
    And grid should have a "cell 1000 of AGE" area
    And grid should not have a "cell 5 of AGE" area
    When user closes the context menu
    Then row 1000 should be current
    And grid should have a "cell 1000 of AGE" area
    And grid should not have a "cell 5 of AGE" area
    And no errors should have been logged

  Scenario: Copy as SMILES copies the right-clicked molecule, not the one that was current
    Given user opens smiles dataset
    Then the table should have 1000 rows
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And "canonical_smiles" of the current row should be "FC(F)(F)c1ccc(OC2CCNCC2)cc1"
    When user scrolls the mouse wheel down 200 times over the "cell 5 of canonical_smiles" area of grid
    Then the "row order" reading of grid should contain "1000"
    When user picks "Current Value > Copy as SMILES" from the context menu of the "cell 1000 of canonical_smiles" area of grid
    Then the clipboard should have the text "CC(C)CCOC(=O)c1ccccc1N"
    And row 1000 should be current
    And "canonical_smiles" of the current row should be "CC(C)CCOC(=O)c1ccccc1N"
    And the "row order" reading of grid should contain "1000"
    And no errors should have been logged

  Scenario: Edit Helm... opens the right-clicked peptide, not the one that was current
    Given user opens helm-peptides dataset
    Then the table should have 540 rows
    When user clicks on the "cell 5 of HELM" area of grid
    Then row 5 should be current
    When user scrolls the mouse wheel down 200 times over the "cell 5 of HELM" area of grid
    Then the "row order" reading of grid should contain "540"
    When user picks "Current Value > Edit Helm..." from the context menu of the "cell 540 of HELM" area of grid
    Then HELM notation tab should be visible
    When user clicks on HELM notation tab
    Then HELM notation should contain the text "PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$$"
    And row 540 should be current
    When user presses Escape
    Then HELM notation tab should be hidden
    And no errors should have been logged
