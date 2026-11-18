@guide @demo
Feature: 2 · Every sequence is a molecule · Sequence structure
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). A real peptide drug, octreotide (Sandostatin) from the BILN_complex sample —
  D-Phe and D-Trp, a Cys2–Cys7 disulfide, a C-terminal threoninol — clicked with the context panel open: Bio's
  Molecular Structure pane draws it atom by atom (Chem's macromolecule assembler); a click on the drawing opens it
  full screen at twice the resolution, the magnifying lens over its disulfide and its indole; the Molecular 3D
  Structure pane shows it in 3D, as large as the widened panel, turned with the pointer; back in the Molecular
  Structure pane, its Explore section holds every Chem panel for the molecule, and the clip ends on its
  Properties (MW 1019, TPSA 383, the formula). The 3D coordinates come from Chem's embed3d, which follows Chem's
  engine setting: Crux on this stand, in the browser (in RDKit mode it runs a Python script, which a demo
  avoids). Set up off camera: the context panel widened (put back at feature end), the columns, the rows and
  monomer names up to eight letters (Thr-ol). Left out: the Composition analysis pane (a table of monomer counts).

  Scenario: A peptide drug in 2D and 3D
    Octreotide's molecule built from its monomers, turned in 3D, explored with Chem's panels.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens BILN_complex dataset keeping the first 16 rows
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is open
    And the context panel is 260 px wider
    And the grid shows only the columns "name, biln"
    And the "biln" column of the grid is 420 px wide
    And the monomers of the "biln" column are drawn at 16 px
    And the monomer names of the "biln" column are cut at 8 characters
    And the rows of the grid are 64 px tall
    When user clicks on the "cell 2 of biln" area of grid
    Then the context panel should show the current cell
    When user expands "Molecular Structure" section
    Then the "Molecular Structure" section should show a drawn molecule
    When the video says "Octreotide, built atom by atom from its monomers"
    And user clicks on the molecule in the "Molecular Structure" section
    And user hovers over the molecule in the "Molecule" dialog at 36% across and 45% down
    And user hovers over the molecule in the "Molecule" dialog at 62% across and 22% down
    And user closes "Molecule" dialog
    And user collapses "Molecular Structure" section
    And user expands "Molecular 3D Structure" section
    Then the "Molecular 3D Structure" section should show a drawn molecule
    When the video says "3D coordinates computed in the browser"
    And user turns the 3D structure in the "Molecular 3D Structure" section
    And user collapses "Molecular 3D Structure" section
    And user expands "Molecular Structure" section
    Then the "Molecular Structure" section should show a drawn molecule
    When the video says "Every chemistry panel works on the peptide"
    And user expands the "Explore" pane inside the "Molecular Structure" section
    And user expands the "Chemistry" pane inside the "Explore" section
    And user expands the "Properties" pane inside the "Chemistry" section
    Then no error or warning balloon should have been shown
