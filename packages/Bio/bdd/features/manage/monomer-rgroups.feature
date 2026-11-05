@journey @serial @realizes:bio.menu.manage.monomers
Feature: The R groups of a monomer in the Manage Monomers editor
  The editor has a chip per R group the drawn monomer carries, under the sketcher, showing its cap; a click
  opens its menu: the fifteen common cap groups, each with what it is and the SMILES Save writes, the current
  one checked, and Edit definition, which opens the fields a library stores under the chips. A library
  definition is read however it is written — [H][*:1], an atom-mapped [O:2], a name where the SMILES belongs
  ([*:3]Azide) — and shown as the cap it means; a cap that is none of the common ones (X) is appended to
  its menu, still flagged while it is a placeholder rather than a cap; an R group nothing defines checks
  nothing. Picking a cap writes it in the libraries' canonical form, and the fields edit the definition as
  typed. While an R group cannot be saved, its chip is red with the reason in its tooltip, and Save is
  disabled — the line beside it saying why. Saving over a monomer the library has asks first.

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
    Then "R1" R-group chip in monomer form should contain the text "R1 H"
    And "R2" R-group chip in monomer form should contain the text "R2 OH"
    And "R3" R-group chip in monomer form should contain the text "R3 Azide"
    And "R4" R-group chip in monomer form should contain the text "R4 X"
    And "R4" R-group chip in monomer form should be invalid
    And "Save Monomer" button in monomer form should be disabled
    And monomer form should contain the text "R4: X is a placeholder, not a cap group — pick the cap R4 gets (usually H)"
    When user clicks on "R4" R-group chip in monomer form
    Then menu should list the item "X"
    And "X" menu item should be checked
    When user presses Escape
    And user hovers over "R5" R-group chip in monomer form
    Then tooltip should contain text "R5: pick its cap group, or edit its definition"
    When user clicks on "R2" R-group chip in monomer form
    Then "OH" menu item should be checked
    And menu should contain text "hydroxyl · O[*:2]"
    When user presses Escape
    Then no errors should have been logged

  Scenario: A pick and the pencil define the R groups, and Save writes them into the library
    When user clicks on "R4" R-group chip in monomer form
    And user clicks on "Br" menu item
    Then "R4" R-group chip in monomer form should be valid
    And "R4" R-group chip in monomer form should contain the text "R4 Br"
    When user clicks on "R5" R-group chip in monomer form
    And user clicks on "Edit definition…" menu item
    And user enters "[*:5]C(F)(F)F" into "Cap SMILES" input in monomer form
    And user enters "CF3" into "Cap name" input in monomer form
    And user enters "R5-CF3" into "Alternate ID" input in monomer form
    Then "R5" R-group chip in monomer form should be valid
    And "R5" R-group chip in monomer form should contain the text "R5 CF3"
    When user clicks on "R5" R-group chip in monomer form
    Then "CF3" menu item should be checked
    When user presses Escape
    Then "Save Monomer" button in monomer form should be enabled
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
