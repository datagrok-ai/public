@guide
Feature: MolTrack EMS: fix what the compound form names
  A guide: the answer to "why will my compound not save, and what do I fix?". Save stays available while
  the form is filled in. Pressed on an incomplete or wrong form, it writes nothing and names every field
  to fix at once, under the field and in the status bar, a required field never touched included. Fix
  them, draw the structure, and the same Save registers the compound.

  Scenario: Fix a mistyped CAS and a missing structure, then register the compound
    # caption: Sign in
    Given user is logged in
    # caption: Simple mode is off
    And simple mode is off
    # caption: Make sure the compound is not registered yet
    And no compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" is registered, now and at feature end
    # caption: Open the Compounds view of MolTrack
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Click New to start a compound
    When user clicks on "New" button in app ribbon
    # caption: Type its common name
    And user types "Triacontanoic acid" into "Common Name" text input
    # caption: Now a mistyped CAS number
    And user types "506-50" into "CAS" text input
    # caption: Press Save
    And user clicks on "Save" button in app ribbon
    # caption: Nothing is saved, and the status bar says why
    Then save problem should contain text "CAS must look like 50-78-2"
    # caption: The CAS field shows what is wrong
    And error of "CAS" text input should be visible
    # caption: The empty Structure is flagged too, though it was never touched
    And error of "Structure" dart input should be visible
    # caption: Nothing reached the server
    And 0 compounds named "Triacontanoic acid" should be registered on the server
    # caption: Fix the CAS number
    When user types "506-50-3" into "CAS" text input
    # caption: Its message goes away
    Then error of "CAS" text input should be hidden
    # caption: Click the structure box to open the sketcher
    When user clicks on molecule canvas in "Structure" dart input
    # caption: Paste the SMILES
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept the structure
    And user clicks on OK button in sketcher dialog
    # caption: The formula is worked out right away
    Then compound identity should contain text "C30H60O2"
    # caption: The Structure is no longer flagged
    And error of "Structure" dart input should be hidden
    # caption: Press Save again
    When user clicks on "Save" button in app ribbon
    # caption: The compound is saved
    Then notification should contain text "Compound saved"
    # caption: It now has its Grok ID
    And breadcrumbs should not contain text "New compound"
    # caption: And it is on the server
    And 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
