@realizes:filters.cp.chem-and-bio-filters @serial
Feature: Substructure filter regressions from closed issues
  Closed issues of the substructure filter card on spgi-100 (100 rows, molecule columns Structure, Core
  and the R groups), each a scenario of its own on a fresh table. With pyridine on the Structure card
  (17 rows), the panel's master switch off lets every row through and on again brings the 17 back
  (public#3328). Current Value | Use as filter on a Core cell filters the Core card (19 rows) and leaves
  the Structure card empty (public#1473, public#1225, public#1216). A search type changed on the
  original view's card is followed by a cloned view, and a structure changed on the clone's card is
  followed by the original (public#2512, public#2639). Use as filter on the first Structure cell still
  filters, to that one molecule, after a cloned view dropped its Structure card (public#2628,
  public#2473). Cancel in the sketcher gives back what the card held before it opened: nothing for a new
  structure, pyridine for an edit (public#697). Pyridine typed into the Structure header's column popup
  is still in its sketcher when the popup is opened again (public#1239, reported on Marvin; here on
  OpenChemLib, the sketcher the feature pins).
  Serial (`@serial`): the Cancel scenarios check Filter as you draw, the account's `chem-filter` setting,
  which is put back at feature end, as card-interplay.feature does.

  Background:
    Given user is logged in
    And the "chem-filter" user settings are put back at feature end
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens spgi-100 dataset
    And user clicks on filter icon in toolbar

  Scenario: The panel switched off and on again filters by the structure the card holds (public#3328)
    When user clicks on "Sketch" text in "Structure" filter card
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user hovers over filter panel
    And user unchecks master of filter panel
    Then the "active" reading of filter panel should be "false"
    And all rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    When user hovers over filter panel
    And user checks master of filter panel
    Then the "active" reading of filter panel should be "true"
    And 17 rows should pass the filter
    And the "filtering of Structure" reading of filter panel should be "true"
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccncc1"
    And no errors should have been logged

  Scenario: Use as filter on a Core cell filters the Core card, not the Structure card (public#1473)
    When user scrolls the grid to the "Core" column
    And user picks "Current Value > Use as filter" from the context menu of the "cell 1 of Core" area of grid
    Then the "filtering of Core" reading of filter panel should be "true"
    And the "structure of Core" reading of filter panel should be the molecule in row 1 of "Core" column
    And the "filtering of Structure" reading of filter panel should be "false"
    And the "structure of Structure" reading of filter panel should be ""
    And 19 rows should pass the filter
    And no errors should have been logged

  Scenario: A search type changed on the original's card is followed by the clone, and a structure changed on the clone by the original (public#2512, public#2639)
    When user clicks on "Sketch" text in "Structure" filter card
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user picks "View > Layout > Clone View" from the top menu
    Then filter panel should be visible
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    When user switches to the "spgi-100" view
    And user opens the settings of the "Structure" filter card
    And user picks search type "Not contains" in the "Structure" filter card
    Then 83 rows should pass the filter
    When user switches to the "spgi-100 copy" view
    Then the "search type of Structure" reading of filter panel should be "Not contains"
    And 83 rows should pass the filter
    When user opens the settings of the "Structure" filter card
    And user picks search type "Contains" in the "Structure" filter card
    Then 17 rows should pass the filter
    When user clicks on the "card Structure" area of filter panel
    And user clears molecule input of sketcher dialog
    And user types "C1CCNCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 15 rows should pass the filter
    When user switches to the "spgi-100" view
    Then the "structure of Structure" reading of filter panel should be "C1CCNCC1"
    And the "search type of Structure" reading of filter panel should be "Contains"
    And 15 rows should pass the filter
    And no errors should have been logged

  Scenario: Use as filter still filters after a cloned view dropped its Structure card (public#2628, public#2473)
    When user picks "View > Layout > Clone View" from the top menu
    Then filter panel should be visible
    When user hovers over "Structure" filter card
    And user clicks on close of "Structure" filter card
    Then "Structure" filter card should be absent
    When user switches to the "spgi-100" view
    And user picks "Current Value > Use as filter" from the context menu of the "cell 1 of Structure" area of grid
    Then the "filtering of Structure" reading of filter panel should be "true"
    And the "structure of Structure" reading of filter panel should be the molecule in row 1 of "Structure" column
    And 1 row should pass the filter
    And no errors should have been logged

  Scenario: Cancel on a new structure drawn with Filter as you draw gives every row back (public#697)
    When user clicks on "Sketch" text in "Structure" filter card
    And user checks "Filter as you draw" input in sketcher dialog
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then 17 rows should pass the filter
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be hidden
    And all rows should pass the filter
    And the "structure of Structure" reading of filter panel should be ""
    And the "filtering of Structure" reading of filter panel should be "false"
    And no errors should have been logged

  Scenario: Cancel on an edited structure drawn with Filter as you draw gives the old structure back (public#697)
    When user clicks on "Sketch" text in "Structure" filter card
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user clicks on the "card Structure" area of filter panel
    And user checks "Filter as you draw" input in sketcher dialog
    And user clears molecule input of sketcher dialog
    And user types "C1CCNCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then 15 rows should pass the filter
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be hidden
    And 17 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And no errors should have been logged

  Scenario: The header popup's substructure filter still holds its molecule when the popup is opened again (public#1239)
    When user hovers over the "header Structure" area of grid
    And user clicks on "Column options" icon in grid
    Then column popup should be visible
    When user types "c1ccncc1" into molecule input of column popup
    And user presses Enter in molecule input of column popup
    Then 17 rows should pass the filter
    When user presses Escape
    Then column popup should be absent
    When user hovers over the "header Structure" area of grid
    And user clicks on "Column options" icon in grid
    Then column popup should be visible
    And the sketcher in column popup should hold the molecule "c1ccncc1"
    And 17 rows should pass the filter
    And no errors should have been logged
