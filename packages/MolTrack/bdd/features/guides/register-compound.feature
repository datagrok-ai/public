@guide
Feature: MolTrack EMS: register a compound
  A guide: the answer to "how do I register a new compound, and what happens to a salt of one that is
  registered already?". Compounds > New opens the compound form; the structure drawn there is
  standardized in the browser (Crux: cleanup, largest fragment, neutralized) and identified by its
  registration hash under the matching rule, so the formula and the identity show under the form
  before Save. Save gives the compound its Grok ID; the sodium salt of the same acid has the same parent,
  so the page names the compound it already is and Save is refused.

  Scenario: Register an acid, then see its sodium salt refused as the same compound
    Given user is logged in
    And simple mode is off
    And no compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" is registered, now and at feature end
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Start a new compound
    When user clicks on "New compound" button in app ribbon
    # caption: Draw the structure: click the empty picture to open the sketcher
    And user clicks on molecule canvas in "Structure" dart input
    # caption: Type or paste a SMILES
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    # caption: The formula and the identity are computed before Save
    Then compound identity should contain text "C29H58O2"
    When user types "Nonacosanoic acid" into "Common Name" text input
    And user clicks on "Save" button in app ribbon
    # caption: Saved: the compound has its Grok ID
    Then notification should contain text "Compound saved"
    And breadcrumbs should not contain text "New compound"
    And 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
    When user remembers the compound the breadcrumbs name
    # caption: Now its sodium salt
    And user clicks on "New compound" button in app ribbon
    And user clicks on molecule canvas in "Structure" dart input
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    # caption: Refused: the salt is the compound registered a moment ago
    Then compound identity should contain text "Already registered as DG-"
    And compound identity should name the remembered compound
    When user clicks on "Save" button in app ribbon
    # caption: Save is refused, and says why
    Then save problem should contain text "Already registered as DG-"
    When user clicks on "Discard" button in app ribbon
    Then 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
