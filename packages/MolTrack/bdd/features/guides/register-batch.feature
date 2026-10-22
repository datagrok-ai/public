@guide
Feature: MolTrack EMS: register a batch of a registered compound
  A guide: the answer to "how do I register a batch?". A batch is a lot of a registered compound. Open the
  compound: its Batches tab lists its batches, and New there starts one with the compound filled in; Save stores
  it. Batches > New is the other way in, picking the compound by its Grok ID number. Many batches at once come
  in through Batches > ⋯ > Import, a column of Grok IDs naming their compounds.

  Scenario: Open the compound, add a batch to it, save it
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And no batch with the EPA Batch ID "GUIDE-001" is registered, now and at feature end
    # caption: Open the compound by its Grok ID: DG-000482, bisphenol A
    And user opens the "Compounds?corporate_compound_id=DG-000482" view of MolTrack EMS
    Then breadcrumbs should contain text "DG-000482"
    # caption: Its batches are below it; New starts one with the compound filled in
    When user clicks on "New" button in batches tab
    # caption: The batch's own fields
    And user types "GUIDE-001" into "EPA Batch ID" text input
    And user types "98" into "Purity (%)" number input
    And user clicks on "Save" button in app ribbon
    # caption: The batch is saved on that compound
    Then 1 batch with the EPA Batch ID "GUIDE-001" should be on the server
