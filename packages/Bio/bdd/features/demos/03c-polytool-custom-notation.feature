@guide @demo
Feature: 1 · Any notation · Custom notation with PolyTool
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). SequenceTranslator's cyclized sample writes cyclic peptides in a
  custom notation, ring closures marked C(1)…C(1); Bio > PolyTool > Convert applies the notation's rules
  (rules_example.json) and writes each sequence as HELM, drawn as its ring, and as a molecule. The pointer on
  a monomer of the custom notation (the `seqs` column: the link to the molecule starts there) lights its atoms
  in the molecule of its row, drawn large in 300 px rows. The rules also run reactions: the clip ends on the
  two click-chemistry rows (7-8), an azide and an alkyne closed into a triazole, in rows taller still, under a
  note and not hovered (a monomer a reaction consumed lights no atoms; most of the reacting monomers' tooltips
  carry polytool-lib's misspelled names, Acetylenil, Ethilenil). Left out: the Diels-Alder and [2+2] rows
  (11-14). The grid is scrolled to rows 7-8 off camera, never filtered (the link looks a row's atoms up by its
  grid row), and the pointer leaves the grid before: a monomer hover the grid still holds when it scrolls is
  drawn on whichever molecule then lies under the pointer. Set up off camera: the columns' widths, the custom
  monomers' names whole (azG(4), C(1)), the row heights and no current row (its color would tint the molecules).

  Scenario: Custom notation to HELM and molecules
    PolyTool reads a custom notation by its rules, from ring closures to click chemistry.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens cyclized dataset
    Then "seqs" column should have semantic type "Macromolecule"
    Given the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the "seqs" column of the grid is 520 px wide
    And the monomers of the "seqs" column are drawn at 15 px
    And the monomer names of the "seqs" column are cut at 8 characters
    And the rows of the grid are 150 px tall
    When user picks "Bio > PolyTool > Convert..." from the top menu
    Then "PolyTool Conversion" dialog should be visible
    When user clicks on OK button in "PolyTool Conversion" dialog
    Then the "PolyTool Conversion" dialog should close
    And a new column "transformed(seqs)" should have been added
    And a new column "molfile(seqs)" should have been added
    Given the grid shows only the columns "seqs, transformed(seqs), molfile(seqs)"
    And the "seqs" column of the grid is 440 px wide
    And the "transformed(seqs)" column of the grid is 320 px wide
    And the "molfile(seqs)" column of the grid is 720 px wide
    And the rows of the grid are 300 px tall
    And no row of the table is current
    When the video says "C(1)…C(1) closures become a HELM ring and macrocycle"
    And user hovers over monomer 3 of the "cell 1 of seqs" area of grid
    And user hovers over monomer 1 of the "cell 2 of seqs" area of grid
    And user moves the pointer away from grid
    And the video says "Click chemistry: azide plus alkyne closes a triazole"
    Given the rows of the grid are 370 px tall
    And the grid is scrolled down to row 7
    Then no errors should have been logged
