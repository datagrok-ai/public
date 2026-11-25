@guide @demo
Feature: 2 · Every sequence is a molecule · Monomers as molecules
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Bio > Transform > Split to Monomers turns the head-to-tail cyclic peptides of
  the HELM_cyclic sample into one column per position, every cell a monomer drawn in its color. A click on a
  non-natural one, 3,4-dichlorophenylalanine (Phe_34diCl), with the context panel open shows Bio's Monomer pane:
  its library record and structure, then — the record folded away — the monomer as a molecule (capped, R-groups
  replaced) with every Chem panel for it, its Properties opened. Set up off camera: the first eight positions
  shown, wide enough for whole names, and the context panel opened right after the click, so it never shows empty.

  Scenario: Monomers are molecules too
    Split into monomers, each non-natural monomer opens with every chemistry panel.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the Bio package is initialized
    And user opens HELM_cyclic dataset keeping the first 8 rows
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence_id, sequence"
    And the "sequence" column of the grid is 260 px wide
    And the rows of the grid are 110 px tall
    When user picks "Bio > Transform > Split to Monomers..." from the top menu
    Then "Split to Monomers" dialog should be visible
    When user clicks on OK button in "Split to Monomers" dialog
    Then the top menu command should have completed
    And 16 new columns should have been added
    And "7" column should have semantic type "Monomer"
    And the value of "7" column in row 1 should be "Phe_34diCl"
    Given the grid shows only the columns "sequence_id, 1, 2, 3, 4, 5, 6, 7, 8"
    And the "1" column of the grid is 96 px wide
    And the "2" column of the grid is 96 px wide
    And the "3" column of the grid is 96 px wide
    And the "4" column of the grid is 96 px wide
    And the "5" column of the grid is 96 px wide
    And the "6" column of the grid is 96 px wide
    And the "7" column of the grid is 96 px wide
    And the "8" column of the grid is 96 px wide
    And the rows of the grid are 36 px tall
    When the video says "One column per position, a monomer per cell"
    And user clicks on the "cell 1 of 7" area of grid
    Given the context panel is open
    Then the context panel should show the current cell
    When user expands "Monomer" section
    And the video says "A monomer is a molecule: every chemistry panel applies"
    And user collapses the "Details" pane inside the "Monomer" section
    And user expands the "Chemistry" pane inside the "Molecule" section
    And user expands the "Properties" pane inside the "Chemistry" section
    Then no error or warning balloon should have been shown
