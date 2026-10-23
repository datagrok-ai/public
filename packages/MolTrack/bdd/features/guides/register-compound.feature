@guide
Feature: MolTrack EMS: register a compound
  A guide: the answer to "how do I register a new compound, and what happens to a salt of one that is
  registered already?". Compounds > New opens the compound form; the structure drawn there is
  standardized in the browser (Crux: cleanup, largest fragment, neutralized) and identified by its
  registration hash under the matching rule, so the formula and the identity show under the form
  before Save. Save gives the compound its Grok ID; the sodium salt of the same acid has the same parent,
  so the Structure field names the compound it already is and Save is refused.

  Scenario: Register an acid, then see its sodium salt refused as the same compound
    # caption: Sign in
    Given user is logged in
    # caption: Simple mode is off
    And simple mode is off
    # caption: Make sure the compound is not registered yet
    And no compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" is registered, now and at feature end
    # caption: Open the Compounds view of MolTrack
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Click New to start a compound
    When user clicks on "New" button in app ribbon
    # caption: Click the structure box to open the sketcher
    And user clicks on molecule canvas in "Structure" dart input
    # caption: Paste the SMILES
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept the structure
    And user clicks on OK button in sketcher dialog
    # caption: The formula and the identity are worked out before Save
    Then compound identity should contain text "C29H58O2"
    # caption: Give it a common name
    When user types "Nonacosanoic acid" into "Common Name" text input
    # caption: Press Save
    And user clicks on "Save" button in app ribbon
    # caption: The compound is saved
    Then notification should contain text "Compound saved"
    # caption: It now has its Grok ID
    And breadcrumbs should not contain text "New compound"
    # caption: And it is on the server
    And 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
    # caption: Remember its Grok ID
    When user remembers the compound the breadcrumbs name
    # caption: Now try its sodium salt: click New
    And user clicks on "New" button in app ribbon
    # caption: Open the sketcher
    And user clicks on molecule canvas in "Structure" dart input
    # caption: Paste the salt
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept the structure
    And user clicks on OK button in sketcher dialog
    # caption: The Structure field says it is the compound registered a moment ago
    Then error of "Structure" dart input should name the remembered compound
    # caption: The identity panel says the same
    And compound identity should contain text "Already registered as DG-"
    # caption: Press Save anyway
    When user clicks on "Save" button in app ribbon
    # caption: Save is refused, and the status bar says why
    Then save problem should contain text "Already registered as DG-"
    # caption: Discard the draft
    When user clicks on "Discard" button in app ribbon
    # caption: Still only one such compound on the server
    Then 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
