@guide
Feature: MolTrack EMS: register compounds in bulk
  A guide: the answer to "how do I register a table of structures?". Compounds > ⋯ > Import opens the import
  wizard: the table (any open table or a file), its columns mapped to the compound's, a preview of every row
  (the structure drawn, and whether it is new, already registered, repeated or unreadable, and why) and a
  report with each row's Grok ID, written back to the table. Batches of these compounds come in through
  Batches > ⋯ > Import, a column of the Grok IDs naming each batch's compound.

  Scenario: Pick a table, see what is new, register it
    # caption: Sign in
    Given user is logged in
    # caption: Simple mode is off
    And simple mode is off
    # caption: The demo data is loaded
    And the MolTrack EMS demo data is on the server
    # caption: Open the example table: three new structures, one already registered, one broken
    And the bulk registration example is open, and its compounds are not registered now and at feature end
    # caption: Open the Compounds view of MolTrack
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Open the more actions menu
    When user clicks on more actions in app ribbon
    # caption: Choose Import
    And user clicks on "Import…" menu item
    # caption: The import wizard opens
    Then "Import compounds" dialog should be visible
    # caption: Pick the table to register
    When user selects "bulk-registration" in "Data" input in "Import compounds" dialog
    # caption: Turn off all or nothing, so one broken row does not stop the rest
    And user unchecks "All or nothing" input in "Import compounds" dialog
    # caption: Next
    And user clicks on "Next" button in "Import compounds" dialog
    # caption: The smiles column fills the structure
    Then import mapping should contain text "smiles"
    # caption: No hash or matching rule to map: they are worked out
    And import mapping should not contain text "Registration Hash"
    # caption: Next
    When user clicks on "Next" button in "Import compounds" dialog
    # caption: The preview: three will be added, one skipped, one with an error
    Then import preview should contain text "3 will be added, 0 updated, 1 skipped, 1 row has errors"
    # caption: Row four is already registered, with its Grok ID
    And import preview should contain text "Already registered as DG-"
    # caption: Row five cannot be read
    And import preview should contain text "Cannot read the structure"
    # caption: Import the three new compounds
    When user clicks on "Import 3 compounds" button in "Import compounds" dialog
    # caption: The report: three inserted, one skipped, one failed
    Then import report should contain text "3 inserted, 0 updated, 1 skipped, 1 failed"
    # caption: Every row has its Grok ID or the reason it has none
    And import report should contain text "Grok ID"
    # caption: Close the wizard
    When user clicks on "Close" button in "Import compounds" dialog
    # caption: The Grok IDs are written back to the table
    Then "Grok ID" column of table "bulk-registration" should have type "string"
    # caption: The three new compounds are on the server
    And the bulk registration example should have 3 new compounds on the server
