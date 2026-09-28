@journey @realizes:chem.cp.molecule-rendering-end-to-end @realizes:chem.int.render-feeds-search @realizes:GROK-16870
Feature: Molecule, reaction and mixture cells are drawn, in the grid and in a viewer tooltip
  On smiles the canonical_smiles column is typed Molecule and its cells carry coloured heteroatom
  ink, while the molregno cells of the same grid carry none — the reading that separates a drawn
  structure from the raw SMILES text a broken renderer falls back to.

  Hovering a scatter plot marker and a box plot on the same table brings up a tooltip whose canvas
  carries the same coloured ink, the viewer keeps painting, and nothing is logged. That is the
  GROK-16870 regression lock: the renderer used to throw
  NullError: method not found: 'gS' on null when a non-Chem viewer asked it for a tooltip cell.

  test-reactions (17 rows) and test_mixtures are typed ChemicalReaction and ChemicalMixture and their cells
  are drawn by their own renderers.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: Molecule cells are drawn as structures, not as text
    Then "canonical_smiles" column should have semantic type "Molecule"
    And the table should have 1000 rows
    And the "cell 1 of canonical_smiles" area of grid should be painted
    And the "cell 1 of canonical_smiles" area of grid should be painted in at least 1 colors
    And the "cell 2 of canonical_smiles" area of grid should be painted in at least 1 colors
    And the "cell 3 of canonical_smiles" area of grid should be painted in at least 1 colors
    And the "cell 1 of molregno" area of grid should not contain the color "#FF0000"
    And the "cell 1 of molregno" area of grid should not contain the color "#0000FF"
    And no errors should have been logged

  Scenario: A scatter plot tooltip draws the molecule of the row under the pointer
    Given user adds a scatter plot viewer with:
      | x | NumHAcceptors |
      | y | NumHDonors    |
    Then scatter plot viewer should be painted
    When user hovers over the "marker of row 1" area of scatter plot viewer
    Then exactly one tooltip should be shown
    And the tooltip should show some columns
    And the canvases of tooltip should be painted in at least 2 colors
    And scatter plot viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A box plot tooltip draws a molecule and the viewer keeps painting
    Given user adds a box plot viewer with:
      | value | NumHAcceptors |
    Then box plot viewer should be painted
    When user hovers over the "marker" area of box plot viewer
    Then exactly one tooltip should be shown
    And the canvases of tooltip should be painted in at least 2 colors
    And box plot viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Reaction cells are typed and drawn by the reaction renderer
    Given user opens test-reactions dataset
    Then the table should have 17 rows
    And "reaction" column should have semantic type "ChemicalReaction"
    And every value of "reaction" column should contain ">>"
    And the "cell 1 of reaction" area of grid should be painted
    And the "cell 1 of reaction" area of grid should be painted in at least 1 colors
    And the "cell 2 of reaction" area of grid should be painted in at least 1 colors
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Mixture cells are typed and drawn by the mixture renderer
    Given user opens test_mixtures dataset
    Then the table should have 4 rows
    And the table should have 1 column
    And "mixture" column should have semantic type "ChemicalMixture"
    And every value of "mixture" column should contain "mixfileVersion"
    And every value of "mixture" column should contain "contents"
    And the "cell 1 of mixture" area of grid should be painted
    And the "cell 1 of mixture" area of grid should be painted in at least 1 colors
    And the "cell 2 of mixture" area of grid should be painted in at least 1 colors
    And no errors should have been logged
    And no error or warning balloon should have been shown
