@guide @demo
Feature: 1 · Any notation · BILN connections
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Peptides of the BILN_complex sample, drawn by Bio's sequence renderer with the
  connections BILN writes in its line: first one of each kind — disulfides, a cyclotide, a head-to-tail ring, a
  staple, two-chain insulin, a dendrimer, a three-arm scaffold — as brackets over the monomers and the chains one
  line each; then four conjugates whose parts are separate chains bonded side by side (Fmoc on a dipeptide, a
  quencher on a spacer on a Lys, octa-Arg carrying a PNA through a beta-Ala, a peptide-drug conjugate on a Cys):
  chains joined end to end are one line with a short line in place of the chain break, at the R-groups' heights,
  and a conjugate on a side chain goes on its own line under its Lys or Cys in the taller rows. Hovering a bonded
  monomer marks its partners and lists its bonds in the tooltip; the pointer works from the right end of each bond,
  so the tooltip opens over the empty end of the rows, not over the bracket or the partner it explains; the pointer
  leaves the grid before the conjugate rows come in, and the clip ends on the last conjugate's tooltip, open over
  the empty rows below. The scaffold's hub is hovered through its last Lys: the hub's own tooltip carries
  polytool-lib's misspelled name (Trihexilamine). Set up off camera: the rows of each part filtered, the columns,
  the row heights and monomer names up to eight letters (Bmt_E, NMe24Abz). Not shown: Show Connections off in the
  column's Renderer Settings (renderer-settings.feature claims the box).

  Scenario: Rings, cross-links and conjugates drawn as connections
    Rings, cross-links and conjugates written in one BILN line, drawn and listed on hover.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens BILN_complex dataset
    Then "biln" column should have semantic type "Macromolecule"
    And "biln" column should have units "biln"
    Given the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And user filters rows where "name" is one of "oxytocin, linaclotide, kalata_B1, cyclosporin_like, stapled_helix, insulin_human, MAP4_NANP, trivalent_cRGDfK_THA"
    And the grid shows only the columns "name, category, biln"
    And the "biln" column of the grid is 1040 px wide
    And the monomers of the "biln" column are drawn at 16 px
    And the monomer names of the "biln" column are cut at 8 characters
    And the rows of the grid are 86 px tall
    When the video says "Disulfides, rings and cross-links: brackets over the monomers"
    And user hovers over monomer 6 of the "cell 1 of biln" area of grid
    And user hovers over monomer 29 of the "cell 33 of biln" area of grid
    And user hovers over monomer 11 of the "cell 35 of biln" area of grid
    And user hovers over monomer 14 of the "cell 42 of biln" area of grid
    And the video says "Two chains, one line each, bonded across"
    And user hovers over monomer 40 of the "cell 48 of biln" area of grid
    And user hovers over monomer 15 of the "cell 62 of biln" area of grid
    And user moves the pointer away from grid
    And the video says "Conjugates: separate chains, bonded side by side"
    Given user filters rows where "name" is one of "Fmoc-FF_hydrogelator, quencher_linker_Lys, R8-bAla-PNA, peptide-drug_conjugate_like"
    And the rows of the grid are 96 px tall
    Then 4 rows should pass the filter
    When user hovers over monomer 2 of the "cell 64 of biln" area of grid
    And user hovers over monomer 6 of the "cell 68 of biln" area of grid
    And user hovers over monomer 9 of the "cell 71 of biln" area of grid
    And user hovers over monomer 17 of the "cell 75 of biln" area of grid
    Then no errors should have been logged
