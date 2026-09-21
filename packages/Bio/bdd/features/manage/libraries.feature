@journey @realizes:bio.menu.manage.monomer-libraries @realizes:bio.op.load_monomer_library @realizes:bio.op.save_monomer_library
Feature: Managing monomer libraries
  Bio | Manage | Monomer Libraries lists every library file on the server with a checkbox. A
  toggle rewrites the user's selection and reloads the monomer library, which Bio announces as
  the bio-monomer-lib-loaded event with the sources it loaded; Add uploads a HELM JSON library
  and loads it, Delete unloads and removes it. The feature starts with every library selected
  and leaves it so. With multiple library storages (the monomerDomainDB package's and the
  files), Add asks which one takes the file; a stand with the files alone skips that dialog.

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

  Scenario: Unchecking a library reloads the monomer library without it, checking it back
    Given user listens for "bio-monomer-lib-loaded" custom event
    When user unchecks "HELMCoreLibrary.json" checkbox
    Then the "bio-monomer-lib-loaded" custom event should have fired
    And the monomer library should not be loaded from "HELMCoreLibrary.json"
    And "HELMCoreLibrary.json" checkbox should be unchecked
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
    And no error or warning balloon should have been shown

  Scenario: Delete unloads the library and removes its file after confirmation
    Given user listens for "bio-monomer-lib-loaded" custom event
    When user clicks on "Delete" icon in "bdd-test-lib.json" checkbox
    Then "Warning" dialog should be visible
    And "Warning" dialog should contain text "bdd-test-lib.json"
    When user clicks on OK button in "Warning" dialog
    Then "bdd-test-lib.json" checkbox should become hidden
    And the "bio-monomer-lib-loaded" custom event should have fired
    And the monomer library should not be loaded from "bdd-test-lib.json"
    And "BDD" should not be a known "PEPTIDE" monomer
    And no error or warning balloon should have been shown
