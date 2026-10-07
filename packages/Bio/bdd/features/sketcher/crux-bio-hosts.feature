@crux @sketcher-controls
Feature: Crux in Bio's hosts
  Bio opens a molecule sketcher in two places this feature claims, with Crux, Chem's own sketcher, pinned as the
  session's: the atomic-level demo, whose peptides' molfiles run to several hundred atoms, opened in the molfile cell's
  editor; and the monomer manager's editor, whose R-groups grid is filled from the R labels of the drawn monomer. What
  Crux holds is read through its status (its "atoms" and "smiles" readings, an area for each atom it draws) and drawn
  on its own controls.

  Needs Chem on the stand. The feature drives Crux's own controls (@sketcher-controls): a run that pins another
  sketcher skips it. Nothing is saved: the demo's table is the page's, and the monomer is never saved.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  # HOST-072
  Scenario: A 339-atom peptide of the atomic-level demo opens whole in Crux, and an edit there lands
    When user calls "Bio:demoBioAtomicLevel" function
    Then the table should have 6 rows
    When user double-clicks on the "cell 1 of molfile(HELM)" area of grid
    Then sketcher dialog should be visible
    And the "atoms" reading of crux sketcher widget should be 339
    And crux sketcher widget should have an "atom 338" area
    When user clicks on crux single bond tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "atoms" reading of crux sketcher widget should be 340
    When user clicks on crux undo button
    Then the "atoms" reading of crux sketcher widget should be 339
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be absent

  # RGROUP-012
  Scenario: R1 and R2 drawn in Crux in the monomer manager fill its R-groups grid with R1 and R2
    Given user opens filter_HELM dataset
    And the Bio package is initialized
    When user picks "Bio > Manage > Monomers" from the top menu
    Then the "Manage Monomers" view should be current
    # the editor keeps the sketcher it made first for the page's life, and its options menu ignores a pick of the
    # session's own sketcher (js-api's Sketcher compares the pick with the session's type, not its own): the editor is
    # switched through OpenChemLib to Crux, whichever it holds
    When user clicks on "Options" icon in monomer sketcher
    And user picks "OpenChemLib" from the open menu
    And user clicks on "Options" icon in monomer sketcher
    And user picks "Crux" from the open menu
    And user clicks on "Add New Monomer" icon
    Then the "ready" reading of crux sketcher widget should be "true"
    And the "atoms" reading of crux sketcher widget should be 0
    When user clicks on crux single bond tool
    And user clicks on crux canvas
    And user clicks on the "atom 1" area of crux sketcher widget
    And user clicks on crux R-group tool
    And user clicks on the "atom 0" area of crux sketcher widget
    And user clicks on crux R1 button
    And user clicks on crux R-Group OK button
    And user clicks on the "atom 2" area of crux sketcher widget
    And user clicks on crux R2 button
    And user clicks on crux R-Group OK button
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]C[*:2]"
    When user clicks on "R-groups" tab
    Then the R-groups grid of the monomer form should list "R1, R2"
