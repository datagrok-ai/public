@journey
Feature: MolTrack EMS: registering a compound
  A chemist registers a structure no one has registered through the compound form, sees its salt
  form refused as the same compound, and finds it again by a structure search. Before that, a Save on
  an incomplete form names every field still to fix, including fields never touched. The form derives the
  identity (formula, hashes) in the browser before Save and shows it under the form; the server
  read-back proves the Save; the search is a structure condition of the filter builder. The compound is deleted when the feature
  ends and swept when it starts.

  Background:
    # caption: Sign in
    Given user is logged in
    # caption: Chem is installed
    And the "Chem" package is installed
    # caption: Make sure the test compound is not registered yet
    And no compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" is registered, now and at feature end
    # caption: Open MolTrack
    And user opens the MolTrack EMS app

  Scenario: Save names every field an incomplete compound still needs, and writes nothing
    # caption: Click New to start a compound
    When user clicks on "New" button in app ribbon
    # caption: A new compound page opens
    Then context panel should not contain text "DG-"
    # caption: Type a common name
    When user types "Octacosanoic acid" into "Common Name" text input
    # caption: Type a CAS number with a mistake
    And user types "123" into "CAS" text input
    # caption: Press Save
    And user clicks on "Save" button in app ribbon
    # caption: Nothing is saved, and the status bar says why
    Then save problem should contain text "CAS must look like 50-78-2"
    # caption: The CAS field shows what is wrong
    And error of "CAS" text input should be visible
    # caption: The empty Structure is flagged too
    And error of "Structure" dart input should be visible
    # caption: Its message reads value can't be empty
    And error of "Structure" dart input should contain text "Value can't be empty"
    # caption: Fix the CAS number
    When user types "506-48-9" into "CAS" text input
    # caption: Press Save again
    And user clicks on "Save" button in app ribbon
    # caption: Now the status bar says the structure is required
    Then save problem should contain text "Structure is required"
    # caption: Nothing reached the server
    And 0 compounds named "Octacosanoic acid" should be registered on the server
    # caption: Discard the draft
    When user clicks on "Discard" button in app ribbon
    # caption: Back on the list
    Then domain list should be visible

  Scenario: A new structure is registered through the form
    # caption: Click New again
    When user clicks on "New" button in app ribbon
    # caption: The page is a new compound
    Then breadcrumbs should contain text "New compound"
    # caption: Open the sketcher
    When user clicks on molecule canvas in "Structure" dart input
    # caption: Paste the SMILES
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept the structure
    And user clicks on OK button in sketcher dialog
    # caption: The sketcher is closed
    Then sketcher dialog should be absent
    # caption: The formula is worked out right away
    And compound identity should contain text "C28H56O2"
    # caption: Type its common name
    When user types "Octacosanoic acid" into "Common Name" text input
    # caption: Press Save
    And user clicks on "Save" button in app ribbon
    # caption: The compound is saved
    Then notification should contain text "Compound saved"
    # caption: It now has its Grok ID
    And breadcrumbs should not contain text "New compound"
    # caption: And it is on the server
    And 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
    # caption: Remember its Grok ID
    When user remembers the compound the breadcrumbs name

  Scenario: Its sodium salt is refused as the compound already registered
    # caption: Now its sodium salt: click New
    When user clicks on "New" button in app ribbon
    # caption: Open the sketcher
    And user clicks on molecule canvas in "Structure" dart input
    # caption: Paste the salt
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept the structure
    And user clicks on OK button in sketcher dialog
    # caption: The identity panel says it is already registered
    Then compound identity should contain text "Already registered as DG-"
    # caption: It names the compound registered a moment ago
    And compound identity should name the remembered compound
    # caption: The Structure field shows the message too
    And error of "Structure" dart input should be visible
    # caption: And names that compound
    And error of "Structure" dart input should name the remembered compound
    # caption: Press Save
    When user clicks on "Save" button in app ribbon
    # caption: Save is refused, and says why
    Then save problem should contain text "Already registered as DG-"
    # caption: Discard the draft
    When user clicks on "Discard" button in app ribbon
    # caption: Back on the list
    Then domain list should be visible
    # caption: Still only one such compound on the server
    And 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
    # caption: The server itself refuses a second compound with the same identity
    And the server should refuse another compound with the remembered compound's registration hash

  Scenario: A substructure condition narrows the list to the compound
    # caption: Note how many compounds the list shows
    When user remembers the text of compound count
    # caption: Add a filter condition
    And user clicks on add condition in filters pane
    # caption: Pick the parent structure
    And user selects "Parent Structure" in condition property
    # caption: Choose has substructure
    And user selects "has substructure" in condition operator
    # caption: Open the sketcher for the condition
    And user clicks on molecule canvas in condition value
    # caption: Paste the substructure
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCC" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept it
    And user clicks on OK button in sketcher dialog
    # caption: One compound is left
    Then compound list should have 1 item in total
    # caption: The count says one compound
    And compound count should have text "1 compound"
    # caption: It is the acid
    And compound list should contain text "Octacosanoic acid"
    # caption: With the Grok ID it was given
    And compound list should name the remembered compound
    # caption: Remove the condition
    When user clicks on condition remove
    # caption: All compounds are back
    Then compound count should have the remembered text
