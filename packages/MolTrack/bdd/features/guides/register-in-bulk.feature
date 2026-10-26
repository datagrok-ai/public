@guide
Feature: MolTrack EMS: register compounds in bulk
  A guide: the answer to "how do I register a table of structures?". Compounds > ⋯ > Import opens the import
  wizard: the table (any open table or a file), its columns mapped to the compound's, a preview of every row
  (the structure drawn, and whether it is new, already registered, repeated or unreadable) and a report with
  each row's Grok ID, written back to the table. Batches of these compounds come in through
  Batches > ⋯ > Import, a column of the Grok IDs naming each batch's compound.

  Scenario: Pick a table, see what is new, register it
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And the bulk registration example is open, and its compounds are not registered now and at feature end
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Bulk registration is the Import of the Compounds view
    When user clicks on more actions in app ribbon
    And user clicks on "Import…" menu item
    Then "Import compounds" dialog should be visible
    # caption: The table to register
    When user selects "bulk-registration" in "Data" input in "Import compounds" dialog
    And user unchecks "All or nothing" input in "Import compounds" dialog
    And user clicks on "Next" button in "Import compounds" dialog
    # caption: Each column of the table, and the compound column it fills
    Then import mapping should contain text "smiles"
    When user clicks on "Next" button in "Import compounds" dialog
    # caption: The preview: three new compounds
    Then import preview should contain text "3 to add"
    And import preview should contain text "1 to skip"
    And import preview should contain text "1 error"
    When user clicks on "Import 3 compounds" button in "Import compounds" dialog
    # caption: The report
    Then import report should contain text "3 added"
    When user clicks on "Close" button in "Import compounds" dialog
    # caption: Every row of the table gets its Grok ID back
    Then "Grok ID" column of table "bulk-registration" should have type "string"
    And the bulk registration example should have 3 new compounds on the server
