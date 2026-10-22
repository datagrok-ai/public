@journey
Feature: MolTrack EMS: registering a compound
  A chemist registers a structure no one has registered through the compound form, sees its salt
  form refused as the same compound, and finds it again by a structure search. The form derives the
  identity (formula, hashes) in the browser before Save; the server read-back proves the Save; the
  search is a structure condition of the filter builder. The compound is deleted when the feature
  ends and swept when it starts.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And no compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" is registered, now and at feature end
    And user opens the MolTrack EMS app

  Scenario: A new structure is registered through the form
    When user clicks on "New" button in app ribbon
    Then breadcrumbs should contain text "New compound"
    When user clicks on molecule canvas in "Structure" dart input
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And "Formula" read-only field should contain text "C28H56O2"
    When user types "Octacosanoic acid" into "Common Name" text input
    And user clicks on "Save" button in app ribbon
    Then notification should contain text "Compound saved"
    And breadcrumbs should not contain text "New compound"
    And 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server
    When user remembers the compound the breadcrumbs name

  Scenario: Its sodium salt is refused as the compound already registered
    When user clicks on "New" button in app ribbon
    And user clicks on molecule canvas in "Structure" dart input
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then "Formula" read-only field should contain text "C28H55NaO2"
    When user clicks on "Save" button in app ribbon
    Then notification should contain text "This structure is already registered as DG-"
    And save problem should name the remembered compound
    And breadcrumbs should contain text "New compound"
    When user clicks on "Discard" button in app ribbon
    Then save problem should be hidden
    And domain list should be visible
    And 1 compound with the structure "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" should be registered on the server

  Scenario: A substructure condition narrows the list to the compound
    When user remembers the text of compound count
    And user clicks on "Show the filters" button
    And user clicks on add condition in filters pane
    And user selects "Parent Structure" in condition property
    And user selects "has substructure" in condition operator
    And user clicks on molecule canvas in condition value
    And user types "CCCCCCCCCCCCCCCCCCCCCCCCCCCC" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then compound list should have 1 item in total
    And compound count should have text "1 compound"
    And compound list should contain text "Octacosanoic acid"
    And compound list should name the remembered compound
    When user clicks on condition remove
    Then compound count should have the remembered text
