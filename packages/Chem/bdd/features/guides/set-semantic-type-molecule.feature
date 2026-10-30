Chem@guide @help:govern/catalog
Feature: Draw an undetected SMILES column as molecules
  A guide: the answer to "if I have a smiles column that is not detected as smiles, how do I set its
  semtype to smiles so that molecules are drawn in grid?". The semantic type of SMILES is Molecule
  (SMILES is its notation, kept in the units tag). Clicking the column header shows the column in the
  Context Panel; its Details pane has a Semantic type field, and entering Molecule there sets the
  type and, since Molecule has one cell renderer, switches the grid to drawing structures.

  The table is reagents whose SMILES are shorter than five characters under a name that does not say
  "smiles", which is what the Chem detector leaves as plain text.

  Scenario: Set the semantic type of a SMILES column to Molecule
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens a table "Reagents" with:
      | Reagent | Equivalents |
      | CCO     | 1.0         |
      | CC#N    | 2.0         |
      | CC=O    | 1.5         |
      | OCCO    | 0.5         |
      | CCCl    | 1.2         |
      | NCCN    | 3.0         |
    # caption: The SMILES are shown as plain text
    And the semantic types of the current table are detected
    And the context panel is open
    # caption: Reagent has no semantic type
    Then "Reagent" column should have semantic type ""
    # caption: Click the Reagent column header
    When user clicks on the "header Reagent" area of grid
    Then the context panel should show "Reagent"
    # caption: Expand Details in the Context Panel
    When user expands Details accordion header in context panel
    # caption: Type Molecule into the Semantic type field
    When user types "Molecule" into "Search or enter semantic type" input in Details pane in context panel
    And user presses Enter
    # caption: Reagent is now a Molecule column
    Then "Reagent" column should have semantic type "Molecule"
    And "Reagent" column should have tag "cell.renderer" equal to "Molecule"
    # caption: The grid draws the molecules
    And the "cell 1 of Reagent" area of grid should be painted in at least 1 colors
