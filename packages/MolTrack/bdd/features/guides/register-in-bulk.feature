@guide
Feature: MolTrack EMS: register compounds and batches in bulk
  A guide: the answer to "how do I register a file of structures with their batches?". Batches > ⋯ >
  Import opens the bulk registration: the table (any open table or a file), the columns mapped onto
  the compound and the batch (matched by name, each one changeable), a preview of what each row
  becomes (a new compound, a structure registered already, or a problem), and the per-row report,
  which opens the batches it registered. "On error" is "Reject all" unless changed: nothing is written
  while any row has a problem; here the bad row is skipped instead.

  Scenario: Map a file onto compounds and batches, preview it, register it
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And the bulk registration example is open, and its compounds are not registered now and at feature end
    And user opens the "Batches" view of MolTrack EMS
    # caption: Bulk registration is the Import of the Batches view
    When user clicks on more actions in app ribbon
    And user clicks on "Import…" menu item
    Then "Register in bulk" dialog should be visible
    # caption: The table, what to register, and what a row with a problem does: here, it is skipped
    And "Data" input in "Register in bulk" dialog should contain text "bulk-registration"
    When user selects "Skip the row" in "On error" input in "Register in bulk" dialog
    And user clicks on "Next" button in "Register in bulk" dialog
    # caption: Each column of the file, and what it registers as
    Then register mapping should contain text "Compound · Common Name"
    And register mapping should contain text "Batch · EPA Batch ID"
    When user clicks on "Next" button in "Register in bulk" dialog
    # caption: The preview: each row by its name and batch ID, and what it becomes
    Then register preview should contain text "6 rows: 3 new compounds, 5 batches (1 of a registered compound), 1 with a problem"
    And register preview should contain text "Sodium hentriacontanoate"
    And register preview should contain text "Same compound as row 1"
    And register preview should contain text "Registered as DG-"
    And register preview should contain text "Cannot read the structure"
    And register preview should contain text "The rows with a problem are skipped."
    When user clicks on "Register" button in "Register in bulk" dialog
    # caption: The report, row by row, with the Grok IDs given
    Then register report should contain text "Registered 5 of 6 rows"
    And register report should contain text "BDD-004"
    And register report should contain text "DGB-"
    # caption: Open the batches just registered
    When user clicks on "Show the 5 batches registered" link in register report
    Then batch list should have 5 items in total
    And batch list should contain text "BDD-001"
    And the bulk registration example should have 3 new compounds and 5 batches on the server
