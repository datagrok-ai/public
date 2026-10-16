@journey @serial @realizes:bio.menu.manage.monomers
Feature: Making and removing a monomer in the Manage Monomers view
  Bio | Manage | Monomers shows the monomers of one library at a time, picked in its Library
  input, beside an editor: a monomer drawn in the editor's sketcher (its SMILES typed into the
  sketcher's molecule field, its attachment points as R groups, a chip each under the sketcher, R1 capped
  with H and R2 with OH as on an amino acid) and given
  a symbol and a name is written into that library by Save, and the library is reloaded with it;
  Delete removes the current monomer after a confirmation. The sketcher is Crux: the editor keeps the
  sketcher it was made with for the page's life, and with OpenChemLib the typed SMILES did not reach
  the drawing (2026-10-13, not looked into further). The library is the feature's own copy of the project's
  bdd-test-lib.json (one peptide monomer, BDD), put on the server under its own name and removed at
  feature end, so no shipped library is touched. Serial: the library selection and the loaded
  library are the account's, as for libraries.feature. Translated from TestTrack
  Bio/bio-manage-libraries-crud.md (the Monomers view) and its spec, which opened the view only.

  What a saved monomer holds (its molfile, R groups, the standardisation the save runs) is the monomer
  manager's package test; the feature claims the save and the delete as the library and the view
  show them.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And user opens filter_HELM dataset
    And the Bio package is initialized
    And all monomer libraries are selected
    And the "fixtures/bdd-test-lib.json" file of the project is the "bdd-monomer-editor.json" monomer library on the server
    When user picks "Bio > Manage > Monomers" from the top menu
    Then the "Manage Monomers" view should be current
    And the monomer sketcher of the Manage Monomers view should be ready

  Scenario: The Library input shows the monomers of the library picked
    When user selects "bdd-monomer-editor" in "Library" input in monomer form
    Then the table should have 1 row
    And the value of "Symbol" column in row 1 should be "BDD"

  Scenario: A monomer drawn in the editor and saved joins the library
    When user clicks on "Add New Monomer" icon
    And user types "CC(C)(C)[C@H](N[1*])C([2*])=O" into molecule input of monomer sketcher
    And user presses Enter in molecule input of monomer sketcher
    Then "R1" R-group chip in monomer form should contain the text "R1 H"
    And "R2" R-group chip in monomer form should contain the text "R2 OH"
    When user enters "BDDTLE" into "Symbol" input in monomer form
    And user enters "bdd tert-leucine" into "Name" input in monomer form
    And user clicks on "Add Monomer" button in monomer form
    Then an info balloon containing "Monomer BDDTLE was successfully saved in library bdd-monomer-editor.json" should have been shown
    And "BDDTLE" should be a known "PEPTIDE" monomer
    And the table should have 2 rows
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Delete removes the current monomer from the library after a confirmation
    Given the value of "Symbol" column in row 2 should be "BDDTLE"
    When user clicks on the "cell 2 of Symbol" area of grid
    Then the "current row" reading of grid should be 2
    When user clicks on "Delete" icon
    Then "Remove Monomer" dialog should be visible
    And "Remove Monomer" dialog should contain text "BDDTLE"
    When user clicks on Remove button in "Remove Monomer" dialog
    Then an info balloon containing "Monomer BDDTLE was successfully removed from bdd-monomer-editor.json library" should have been shown
    And "BDDTLE" should not be a known "PEPTIDE" monomer
    And the table should have 1 row
    And no error or warning balloon should have been shown
    And no errors should have been logged
