@journey @serial @realizes:bio.menu.manage.monomers
Feature: The R groups of a monomer in the Manage Monomers editor
  The editor has a chip per R group the drawn monomer carries, under the sketcher, showing its cap; a click
  opens the R group's row above the chips: a choice of the fifteen common cap groups, a line saying what the
  chosen cap is with the SMILES Save writes, and a pencil that opens the four fields a library stores. A library definition is read however it is written —
  [H][*:1], an atom-mapped [O:2], a name where the SMILES belongs ([*:3]Azide) — and shown as the
  cap it means; a cap that is none of the common ones (X) is appended to its row's list, still
  flagged while it is a placeholder rather than a cap; an R group nothing defines shows the empty choice. Picking a
  cap writes it in the libraries' canonical form, and the four fields edit the definition as typed.
  While an R group cannot be saved, its chip is red, its row says why, and Save is disabled — the line
  beside it saying why. Saving over a monomer the library has asks first.

  The library is the feature's own copy of the project's bdd-rgroups-lib.json (one peptide monomer,
  bddRG, with its five R groups defined five ways), put on the server under its own name and removed
  at feature end, so no shipped library is touched. Serial: the library selection and the loaded
  library are the account's. How a definition is read (every form the production libraries use) is
  the package test's (`src/tests/monomer-manager-tests.ts`); the feature claims what the user sees and
  what Save writes.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And user opens filter_HELM dataset
    And the Bio package is initialized
    And all monomer libraries are selected
    And the "fixtures/bdd-rgroups-lib.json" file of the project is the "bdd-rgroups.json" monomer library on the server
    When user picks "Bio > Manage > Monomers" from the top menu
    Then the "Manage Monomers" view should be current
    And the monomer sketcher of the Manage Monomers view should be ready
    When user selects "bdd-rgroups" in "Library" input in monomer form
    Then the table should have 1 row

  Scenario: Each R group shows the cap its definition means, however it is written
    When user clicks on the "cell 1 of Symbol" area of grid
    And user clicks on "Edit Monomer" icon
    Then "R1" input in monomer form should have value "H"
    And "R2" input in monomer form should have value "OH"
    And "R3" input in monomer form should have value "Azide"
    And "R4" R-group chip in monomer form should contain the text "X"
    When user clicks on "R4" R-group chip in monomer form
    Then "R4" input in monomer form should list the item "X"
    And "R4" input in monomer form should have the error "X is a placeholder, not a cap group — pick the cap R4 gets (usually H)"
    And "Save Monomer" button in monomer form should be disabled
    And monomer form should contain the text "R4: X is a placeholder, not a cap group — pick the cap R4 gets (usually H)"
    When user clicks on "R5" R-group chip in monomer form
    Then "R5" input in monomer form should have the error "Pick its cap group, or define it with the pencil"
    When user clicks on "R2" R-group chip in monomer form
    Then "R2" input in monomer form should contain the text "hydroxyl · O[*:2]"
    And no errors should have been logged

  Scenario: A pick and the pencil define the R groups, and Save writes them into the library
    When user clicks on "R4" R-group chip in monomer form
    And user selects "Br" in "R4" input in monomer form
    Then "R4" input in monomer form should have no error
    And "R4" input in monomer form should contain the text "bromo · Br[*:4]"
    And "R4" R-group chip in monomer form should contain the text "Br"
    When user clicks on "R5" R-group chip in monomer form
    And user clicks on "Edit the R5 definition" icon button
    And user enters "[*:5]C(F)(F)F" into "Cap SMILES" input in monomer form
    And user enters "CF3" into "Cap name" input in monomer form
    And user enters "R5-CF3" into "Alternate ID" input in monomer form
    Then "R5" input in monomer form should have no error
    And "R5" input in monomer form should list the item "CF3"
    And "Save Monomer" button in monomer form should be enabled
    When user clicks on "Save Monomer" button in monomer form
    Then "Save Monomer" dialog should be visible
    And "Save Monomer" dialog should contain text "Monomer with symbol 'bddRG' already exists in library bdd-rgroups.json"
    When user clicks on Save button in "Save Monomer" dialog
    And an info balloon containing "Monomer bddRG was successfully saved in library bdd-rgroups.json" should have been shown
    And the value of "R1" column in row 1 should be "[*:1][H]"
    And the value of "R2" column in row 1 should be "O[*:2]"
    And the value of "R3" column in row 1 should be "[N-]=[N+]=N[*:3]"
    And the value of "R4" column in row 1 should be "Br[*:4]"
    And the value of "R5" column in row 1 should be "[*:5]C(F)(F)F"
    And no error or warning balloon should have been shown
    And no errors should have been logged
