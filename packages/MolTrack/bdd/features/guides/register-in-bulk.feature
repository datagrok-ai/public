@guide
Feature: MolTrack EMS: register compounds in bulk
  A guide: the answer to "how do I register a file of structures?". Compounds > ⋯ > Import opens the import
  wizard: the table (any open table or a file), the columns mapped onto the compound (matched by name, a column
  of structures by its type, each one changeable), a preview of what each row becomes (added, skipped as a
  structure registered already, or an error), and the per-row report with each new compound's Grok ID, which
  shows the compounds it registered. "All or nothing" is on unless changed: nothing is written while any row
  has an error; here the bad row fails alone instead. A structure registered already is skipped either way.
  Batches of these compounds come in through Batches > ⋯ > Import, a column of the Grok IDs naming each
  batch's compound.

  Scenario: Map a file onto compounds, preview it, register it
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And the bulk registration example is open, and its compounds are not registered now and at feature end
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Bulk registration is the Import of the Compounds view
    When user clicks on more actions in app ribbon
    And user clicks on "Import…" menu item
    Then "Import compounds" dialog should be visible
    # caption: The table, and what a row with an error does: with "All or nothing" off, it fails alone
    And "Data" input in "Import compounds" dialog should contain text "bulk-registration"
    When user unchecks "All or nothing" input in "Import compounds" dialog
    And user clicks on "Next" button in "Import compounds" dialog
    # caption: Each column of the file, and the compound column it fills
    Then import mapping should contain text "Common Name"
    When user clicks on "Next" button in "Import compounds" dialog
    # caption: The preview: what each row becomes
    Then import preview should contain text "3 will be added, 0 updated, 1 skipped, 1 row has errors."
    And import preview should contain text "This structure is already registered as DG-"
    And import preview should contain text "Cannot read the structure"
    When user clicks on "Next" button in "Import compounds" dialog
    # caption: The report, with each new compound's Grok ID
    Then import report should contain text "3 inserted, 0 updated, 1 skipped, 1 failed."
    And import report should contain text "DG-"
    # caption: Show the compounds just registered
    When user clicks on "Show the 3 imported rows" link in import report
    Then compound list should have 3 items in total
    And compound list should contain text "Hentriacontanoic acid"
    And the bulk registration example should have 3 new compounds on the server
