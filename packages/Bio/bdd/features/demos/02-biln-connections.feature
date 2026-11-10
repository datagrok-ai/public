@guide @demo
Feature: Bioinformatics · BILN connections
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Eight peptides of the BILN_complex sample, one of each kind of
  connection BILN writes in its line — disulfides, a cyclotide, a head-to-tail ring, a staple, two-chain
  insulin, a dendrimer, a three-arm scaffold — drawn by Bio's sequence
  renderer with the connections as brackets over the monomers and the chains one line each. Hovering a
  bonded monomer marks its partners and lists its bonds in the tooltip. Not shown: Show Connections off
  in the column's Renderer Settings (renderer-settings.feature claims the box).

  Scenario: Rings, cross-links and chains drawn as connections
    Cyclizations, disulfides and bonds between chains, written in one BILN line, are drawn as brackets over the monomers and listed on hover.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens BILN_complex dataset
    Then "biln" column should have semantic type "Macromolecule"
    And "biln" column should have units "biln"
    Given user filters rows where "name" is one of "oxytocin, linaclotide, kalata_B1, cyclosporin_like, stapled_helix, insulin_human, MAP4_NANP, trivalent_cRGDfK_THA"
    And the grid shows only the columns "name, category, biln"
    And the "biln" column of the grid is 1040 px wide
    And the monomers of the "biln" column are drawn at 16 px
    And the rows of the grid are 86 px tall
    When user hovers over monomer 1 of the "cell 1 of biln" area of grid
    And the video says "Disulfides, rings and cross-links: brackets over the monomers"
    And user hovers over monomer 5 of the "cell 4 of biln" area of grid
    And user hovers over monomer 1 of the "cell 33 of biln" area of grid
    And user hovers over monomer 4 of the "cell 42 of biln" area of grid
    And the video says "Two chains, one line each, bonded across"
    And user hovers over monomer 7 of the "cell 48 of biln" area of grid
    And user hovers over monomer 9 of the "cell 61 of biln" area of grid
    And user hovers over monomer 16 of the "cell 62 of biln" area of grid
    Then no errors should have been logged
