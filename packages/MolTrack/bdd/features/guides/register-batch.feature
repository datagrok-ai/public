@guide
Feature: MolTrack EMS: register a batch of a registered compound
  A guide: the answer to "how do I register a batch?". A batch is a lot of a registered compound: Batches >
  New opens the batch form, its Compound field picks the compound by its Grok ID or its name, and the
  compound's card shows above the form. A compound that is not registered yet is registered under Compounds
  first; many batches at once come in through Batches > ⋯ > Import, a column of Grok IDs naming their compounds.

  Scenario: Pick the compound, fill in the batch, save it
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And no batch with the EPA Batch ID "GUIDE-001" is registered, now and at feature end
    And user opens the "Batches" view of MolTrack EMS
    # caption: Start a new batch
    When user clicks on "New" button in app ribbon
    # caption: Pick its compound by Grok ID or name
    And user types "Bisphenol A" into "Compound" domain pick
    And user presses Enter in "Compound" domain pick
    # caption: The compound the batch is of
    Then compound card should contain text "Bisphenol A"
    # caption: The batch's own fields
    When user types "GUIDE-001" into "EPA Batch ID" text input
    And user types "98.2" into "Purity" number input
    And user clicks on "Save" button in app ribbon
    # caption: The batch is saved on that compound
    Then notification should contain text "saved"
    And compound card should contain text "Bisphenol A"
    And 1 batch with the EPA Batch ID "GUIDE-001" should be on the server
