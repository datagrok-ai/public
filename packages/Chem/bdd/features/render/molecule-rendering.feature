@journey @realizes:chem.cp.molecule-rendering-end-to-end @realizes:chem.int.render-feeds-search @realizes:GROK-16870
Feature: A molecule drawn in a viewer tooltip, and mixture cells
  Hovering a scatter plot marker and a box plot on smiles brings up a tooltip whose canvas carries
  coloured molecule ink, the viewer keeps painting, and nothing is logged. That is the
  GROK-16870 regression lock: the renderer used to throw
  NullError: method not found: 'gS' on null when a non-Chem viewer asked it for a tooltip cell.

  test_mixtures is typed ChemicalMixture and its cells are drawn by its own renderer.

  Molecule and reaction cells in the grid are tested in Chem src/tests/rendering-tests.ts
  ('rdkit grid cell renderer') and src/tests/reaction-rendering-tests.ts ('Reaction rendering').

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

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
