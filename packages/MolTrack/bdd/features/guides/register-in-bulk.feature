@guide
Feature: MolTrack EMS: register compounds in bulk
  A guide: the answer to "how do I register a file of structures?". Compounds > ⋯ > Import opens the bulk
  registration: the table (any open table or a file), the columns mapped onto the compound (matched by name,
  each one changeable), a preview of what each row becomes (a new compound, a structure registered already, or
  a problem), and the per-row report with each row's Grok ID, which opens the compounds it registered. "On
  error" is "Reject all" unless changed: nothing is written while any row has a problem; here the bad row is
  skipped instead. Batches of these compounds come in through Batches > ⋯ > Import, a column of the Grok IDs
  naming each batch's compound.

  Scenario: Map a file onto compounds, preview it, register it
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And the bulk registration example is open, and its compounds are not registered now and at feature end
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Bulk registration is the Import of the Compounds view
    When user clicks on more actions in app ribbon
    And user clicks on "Import…" menu item
    Then "Register compounds in bulk" dialog should be visible
    # caption: The table, and what a row with a problem does: here, it is skipped
    And "Data" input in "Register compounds in bulk" dialog should contain text "bulk-registration"
    When user selects "Skip the row" in "On error" input in "Register compounds in bulk" dialog
    And user clicks on "Next" button in "Register compounds in bulk" dialog
    # caption: Each column of the file, and what it registers as
    Then register mapping should contain text "Common Name"
    When user clicks on "Next" button in "Register compounds in bulk" dialog
    # caption: The preview: each row by its name, and what it becomes
    Then register preview should contain text "6 rows: 3 new compounds, 1 repeating an earlier row, 1 registered already, 1 with a problem"
    And register preview should contain text "Sodium hentriacontanoate"
    And register preview should contain text "Same compound as row 1"
    And register preview should contain text "This structure is already registered as DG-"
    And register preview should contain text "Cannot read the structure"
    And register preview should contain text "The rows with a problem are skipped."
    When user clicks on "Register" button in "Register compounds in bulk" dialog
    # caption: The report, row by row, with each row's Grok ID
    Then register report should contain text "Registered 3 of 6 rows"
    And register report should contain text "Same compound as row 1"
    And register report should contain text "DG-"
    # caption: Open the compounds just registered
    When user clicks on "Show the 3 compounds registered" link in register report
    Then compound list should have 3 items in total
    And compound list should contain text "Hentriacontanoic acid"
    And the bulk registration example should have 3 new compounds on the server
