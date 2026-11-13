@guide @demo @viewport:1440x810
Feature: Bioinformatics · Sequence structure
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). A head-to-tail cyclic HELM peptide of non-natural monomers, clicked with the
  context panel open: Bio's Molecular Structure pane draws it atom by atom (Chem's macromolecule assembler), a
  click on the drawing opens it full screen at twice the resolution with a magnifying lens under the pointer, and
  the Molecular 3D Structure pane shows it in 3D, turned with the pointer. The 3D coordinates come from Chem's
  embed3d, which follows Chem's engine setting: Crux on this stand, in the browser (in RDKit mode it runs a Python
  script, which a demo avoids). Left out: the Composition analysis pane (a table of monomer counts).

  Scenario: A sequence's structure in 2D and 3D
    Click a peptide: its molecule, built from the monomer library, opens full screen under a magnifier and turns in 3D.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens HELM_cyclic dataset keeping the first 8 rows
    And the toolbox pane is hidden
    And the context panel is open
    And the grid shows only the columns "sequence_id, sequence, activity"
    And the "sequence" column of the grid is 400 px wide
    And the rows of the grid are 190 px tall
    When user clicks on the "cell 2 of sequence" area of grid
    Then the context panel should show the current cell
    When user expands "Molecular Structure" section
    Then the "Molecular Structure" section should show a drawn molecule
    When the video says "Built atom by atom from the monomer library"
    And user clicks on the molecule in the "Molecular Structure" section
    And user hovers over the molecule in the "Molecule" dialog at 32% across and 30% down
    And user hovers over the molecule in the "Molecule" dialog at 62% across and 28% down
    And user hovers over the molecule in the "Molecule" dialog at 58% across and 75% down
    And user closes "Molecule" dialog
    And user expands "Molecular 3D Structure" section
    Then the "Molecular 3D Structure" section should show a drawn molecule
    When the video says "3D coordinates computed in the browser"
    And user turns the 3D structure in the "Molecular 3D Structure" section
    Then no error or warning balloon should have been shown
