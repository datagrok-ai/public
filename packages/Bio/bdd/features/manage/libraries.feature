@journey @serial @realizes:bio.menu.manage.monomer-libraries @realizes:bio.op.load_monomer_library @realizes:bio.op.save_monomer_library
Feature: Managing monomer libraries
  Bio | Manage | Monomer Libraries lists every library file on the server with a checkbox. A
  toggle rewrites the user's selection and reloads the monomer library, which Bio announces as
  the bio-monomer-lib-loaded event with the sources it loaded; Add uploads a HELM JSON library
  and loads it, Delete unloads and removes it. The feature starts with every library selected
  and leaves it so. With multiple library storages (the monomerDomainDB package's and the
  files), Add asks which one takes the file; a stand with the files alone skips that dialog. The
  shipped library file is a HELM library (every monomer has a symbol and a structure); an uploaded
  library is still listed when the manager opens again, and Delete removes the file from every
  storage, not only the checkbox. The dialog entry (Bio:manageMonomerLibraries, the panel link)
  lists the same libraries as the view. Serial: the library selection is the user's, and other
  features (atomic level, the HELM project round-trip) depend on it.

  Not translated: writing a library file straight to the file share and editing it there (the
  lifecycle md's S2.1-2.2 — the old spec wrote a file and read it back, which tests the file
  store, not Bio; the upload through Add is the product path and is claimed); the Manage Monomers
  view's CRUD (creating and editing a monomer through its editor dialog — no feature yet).

  Background:
    Given user is logged in
    And user opens filter_HELM dataset
    And the Bio package is initialized
    And all monomer libraries are selected
    And no "bdd-test-lib.json" monomer library is on the server
    When user picks "Bio > Manage > Monomer Libraries" from the top menu
    Then the top menu command should have completed
    And the "Manage Monomer Libraries" view should be current
    And "Manage Duplicate Monomer Symbols" heading should be visible
    And "HELMCoreLibrary.json" checkbox should be checked
    And Search input should be visible
    And Add button should be visible
    And Merge button should be visible

  Scenario: The shipped library file is a HELM monomer library
    Then the "HELMCoreLibrary.json" monomer library file should list monomers with a symbol and a structure each
    And the monomer library should be loaded from "HELMCoreLibrary.json"

  Scenario: Unchecking a library reloads the monomer library without it, checking it back
    Given user listens for "bio-monomer-lib-loaded" custom event
    When user unchecks "HELMCoreLibrary.json" checkbox
    Then the "bio-monomer-lib-loaded" custom event should have fired
    And the monomer library should not be loaded from "HELMCoreLibrary.json"
    And "HELMCoreLibrary.json" checkbox should be unchecked
    And "polytool-lib.json" checkbox should be checked
    And the monomer library should be loaded from "polytool-lib.json"
    When user checks "HELMCoreLibrary.json" checkbox
    Then the "bio-monomer-lib-loaded" custom event should have fired
    And the monomer library should be loaded from "HELMCoreLibrary.json"
    And "A" should be a known "PEPTIDE" monomer
    And no error or warning balloon should have been shown

  Scenario: Add uploads a library file and its monomers become known
    Given user listens for "bio-monomer-lib-loaded" custom event
    And "BDD" should not be a known "PEPTIDE" monomer
    When user uploads "fixtures/bdd-test-lib.json" through Add button
    And user chooses "Files" storage for the uploaded monomer library
    Then "bdd-test-lib.json" checkbox should become visible
    And "bdd-test-lib.json" checkbox should be checked
    And the "bio-monomer-lib-loaded" custom event should have fired
    And the monomer library should be loaded from "bdd-test-lib.json"
    And "BDD" should be a known "PEPTIDE" monomer
    And an info balloon containing "Added bdd-test-lib.json HELM library" should have been shown
    And no error or warning balloon should have been shown

  Scenario: The uploaded library is still listed when the manager opens again
    When user closes the current view
    And user picks "Bio > Manage > Monomer Libraries" from the top menu
    Then the "Manage Monomer Libraries" view should be current
    And "bdd-test-lib.json" checkbox should be visible
    And "bdd-test-lib.json" checkbox should be checked
    And "HELMCoreLibrary.json" checkbox should be checked

  Scenario: Delete unloads the library and removes its file after confirmation
    Given user listens for "bio-monomer-lib-loaded" custom event
    When user clicks on "Delete" icon in "bdd-test-lib.json" checkbox
    Then "Warning" dialog should be visible
    And "Warning" dialog should contain text "bdd-test-lib.json"
    When user clicks on OK button in "Warning" dialog
    Then "bdd-test-lib.json" checkbox should become absent
    And there should be no "bdd-test-lib.json" monomer library on the server
    And the "bio-monomer-lib-loaded" custom event should have fired
    And the monomer library should not be loaded from "bdd-test-lib.json"
    And "BDD" should not be a known "PEPTIDE" monomer
    And no error or warning balloon should have been shown

  Scenario: The dialog entry lists the same libraries as the view
    When user closes the current view
    And user calls "Bio:manageMonomerLibraries" function
    Then "Manage monomer libraries" dialog should be visible
    And "HELMCoreLibrary.json" checkbox in "Manage monomer libraries" dialog should be checked
    And "polytool-lib.json" checkbox in "Manage monomer libraries" dialog should be checked
    When user presses Escape
    Then "Manage monomer libraries" dialog should be absent

  Scenario: After the dialog, the view lists the libraries again
    When user picks "Bio > Manage > Monomer Libraries" from the top menu
    Then the "Manage Monomer Libraries" view should be current
    And "HELMCoreLibrary.json" checkbox should be visible
    And no error or warning balloon should have been shown
