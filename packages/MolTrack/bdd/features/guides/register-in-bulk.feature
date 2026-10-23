@guide
Feature: MolTrack EMS: register compounds in bulk
  A guide: the answer to "how do I register a table of structures?". Compounds > ⋯ > Import asks for the
  table (any open table or a file), works out each row's identity from its structure and previews it: the
  structure, the verdict (new, already registered, repeated, unreadable) and why. The new compounds go on to
  the import wizard (columns, preview, report), and every row of the table gets its Grok ID back. Batches of these compounds come in through Batches > ⋯ > Import, a column of the
  Grok IDs naming each batch's compound.

  Scenario: Pick a table, see what is new, register it
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And the bulk registration example is open, and its compounds are not registered now and at feature end
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Bulk registration is the Import of the Compounds view
    When user clicks on more actions in app ribbon
    And user clicks on "Import…" menu item
    Then "Register compounds" dialog should be visible
    # caption: The table to register
    And "Data" input in "Register compounds" dialog should contain text "bulk-registration"
    # caption: The preview: each row's structure and what it becomes
    And "Register compounds" dialog should contain text "5 rows: 3 new, 1 already registered, 1 cannot be read"
    When user clicks on "OK" button in "Register compounds" dialog
    # caption: The new compounds go on to the import: columns, preview, report
    Then "Import compounds" dialog should be visible
    When user clicks on "Next" button in "Import compounds" dialog
    # caption: Each column of the table, and the compound column it fills
    Then import mapping should contain text "Common Name"
    When user clicks on "Next" button in "Import compounds" dialog
    # caption: The preview: three new compounds
    Then import preview should contain text "3 will be added"
    When user clicks on "Next" button in "Import compounds" dialog
    # caption: The report
    Then import report should contain text "3 inserted"
    When user clicks on "Close" button in "Import compounds" dialog
    # caption: Every row of the table gets its Grok ID back
    Then "Grok ID" column of table "bulk-registration" should have type "string"
    And the bulk registration example should have 3 new compounds on the server
