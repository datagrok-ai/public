@journey @realizes:bio.app.monomer-collections
Feature: Monomer collections
  The Monomer Collections app shows the collection files under the package's
  monomer-collections folder as cards, the stock ones and the user's. A collection is made from
  the New Collection card and lands on the server as it was entered; a card selects on click;
  Delete removes the file after confirmation.

  Background:
    Given user is logged in
    And the Bio package is initialized
    And no "bdd-test-collection" monomer collection is on the server
    And user opens the "Monomer Collections" app
    Then "Canonical AAs" card should be visible
    And "New Collection" card should be visible

  Scenario: A new collection is made from the dialog and saved on the server
    When user clicks on "New Collection" card
    Then "New Monomer Collection" dialog should be visible
    And "Polymer Type" input in "New Monomer Collection" dialog should have value "PEPTIDE"
    When user types "bdd-test-collection" into Name input in "New Monomer Collection" dialog
    And user types "made by the bdd feature" into Description input in "New Monomer Collection" dialog
    And user types "A" into Search input in "New Monomer Collection" dialog
    And user presses Enter in Search input in "New Monomer Collection" dialog
    And user types "G" into Search input in "New Monomer Collection" dialog
    And user presses Enter in Search input in "New Monomer Collection" dialog
    Then "New Monomer Collection" dialog should contain text "2 monomer(s) selected"
    When user clicks on OK button in "New Monomer Collection" dialog
    Then "New Monomer Collection" dialog should be hidden
    And "bdd-test-collection" card should become visible
    And "bdd-test-collection" card should contain text "2 monomer(s)"
    And "bdd-test-collection" card should contain text "made by the bdd feature"
    And the "bdd-test-collection" monomer collection should hold monomers "A, G"
    And no error or warning balloon should have been shown

  Scenario: A card selects on click
    When user clicks on "bdd-test-collection" card
    Then "bdd-test-collection" card should be selected
    And "Canonical AAs" card should not be selected

  Scenario: Delete removes the collection after confirmation
    When user clicks on Delete button in "bdd-test-collection" card
    Then "Delete Collection" dialog should be visible
    And "Delete Collection" dialog should contain text "bdd-test-collection"
    When user clicks on OK button in "Delete Collection" dialog
    Then "bdd-test-collection" card should become hidden
    And there should be no "bdd-test-collection" monomer collection on the server
    And no error or warning balloon should have been shown
