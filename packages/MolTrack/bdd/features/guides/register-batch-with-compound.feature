@guide
Feature: MolTrack EMS: register a batch together with its compound
  A guide: the answer to "how do I register a batch of a compound that is not registered yet?". Batches >
  New opens the batch form with a "New compound" form above it: drawing a structure there points the
  batch at that new compound, and one Save writes both in one transaction. When the structure is
  registered already, the batch is added to that compound instead, as the legacy registration did.

  Scenario: A new batch registers its new compound in the same save, a second batch joins it
    Given user is logged in
    And simple mode is off
    And no compound with the structure "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO" is registered, now and at feature end
    And user opens the "Batches" view of MolTrack EMS
    # caption: Start a new batch: the form for its compound comes with it
    When user clicks on "New" button in app ribbon
    Then new compound section should be visible
    When user clicks on molecule canvas in new compound section
    And user types "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    And user types "Octacosane-1,28-diol" into "Common Name" text input in new compound section
    # caption: The batch's own fields
    And user types "GUIDE-001" into "EPA Batch ID" text input
    And user types "98.2" into "Purity" number input
    And user clicks on "Save" button in app ribbon
    # caption: One save registered the compound and its batch
    Then notification should contain text "saved"
    And 1 compound with the structure "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO" should be registered on the server
    And the compound with the structure "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO" should have 1 batch on the server
    # caption: A second batch of the same structure
    When user clicks on "New" button in app ribbon
    And user clicks on molecule canvas in new compound section
    And user types "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    # caption: The structure is registered: the batch is added to that compound
    Then registered compound note should contain text "the batch is added to it"
    And new compound section should be hidden
    When user types "GUIDE-002" into "EPA Batch ID" text input
    And user clicks on "Save" button in app ribbon
    Then 1 compound with the structure "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO" should be registered on the server
    And the compound with the structure "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO" should have 2 batches on the server
