@journey @realizes:filters.cp.chem-and-bio-filters
Feature: A substructure card from a cell, after a reset and on a cloned view
  On spgi-100, Current Value | Use as filter on a Structure cell adds one Structure card holding that
  cell's molecule (GROK-20001). A Structure card rebuilt through Remove All and a reopened panel,
  holding pyridine in Similar mode, is cleared by the panel's reset icon: every row passes, the card
  stays, its structure is empty, its search type is Contains, it no longer filters and the header
  counter is hidden, and the card offers its Sketch link again (GROK-14028, GROK-20739). Set again to Similar at cutoff 0.6, the card comes back in a view
  cloned with View | Layout | Clone View as Similar at 0.6 with the same fingerprint and molecule
  (GROK-18530).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens spgi-100 dataset

  Scenario: Use as filter adds one card holding the cell's molecule
    When user clicks on filter icon in toolbar
    And user hovers over "Structure" filter card
    And user clicks on close of "Structure" filter card
    Then "Structure" filter card should be absent
    When user picks "Current Value > Use as filter" from the context menu of the "cell 1 of Structure" area of grid
    Then there should be 1 visible "Structure" filter card
    And the "type of Structure" reading of filter panel should be "Chem:substructureFilter"
    And the "structure of Structure" reading of filter panel should be the molecule of row 1 of "Structure" column
    And fewer than 100 rows should pass the filter
    And the "filtering of Structure" reading of filter panel should be "true"
    And no errors should have been logged

  Scenario: The reset icon clears the card's structure and search type and keeps the card
    When user picks "Remove All" from the viewer menu of filter panel
    And user clicks on close icon of filters viewer
    And user clicks on filter icon in toolbar
    Then there should be 1 visible "Structure" filter card
    And the "structure of Structure" reading of filter panel should be ""
    When user clicks on "Sketch" text in "Structure" filter card
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user opens the settings of the "Structure" filter card
    And user picks search type "Similar" in the "Structure" filter card
    Then the "search type of Structure" reading of filter panel should be "Similar"
    And counter of filter panel should have text "1"
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And "Structure" filter card should be visible
    And the "type of Structure" reading of filter panel should be "Chem:substructureFilter"
    And the "structure of Structure" reading of filter panel should be ""
    And the "search type of Structure" reading of filter panel should be "Contains"
    And the "filtering of Structure" reading of filter panel should be "false"
    And counter of filter panel should be hidden
    And "Sketch" text in "Structure" filter card should be visible
    And no errors should have been logged

  Scenario: A Similar card at a non-default cutoff comes back the same in a cloned view
    When user clicks on "Sketch" text in "Structure" filter card
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user opens the settings of the "Structure" filter card
    And user picks search type "Similar" in the "Structure" filter card
    And user sets the similarity cutoff of the "Structure" filter card to 0.6
    Then the "search type of Structure" reading of filter panel should be "Similar"
    And the "similarity cutoff of Structure" reading of filter panel should be 0.6
    And the "fingerprint of Structure" reading of filter panel should be "Morgan"
    When user picks "View > Layout > Clone View" from the top menu
    Then filter panel should be visible
    And the "search type of Structure" reading of filter panel should be "Similar"
    And the "similarity cutoff of Structure" reading of filter panel should be 0.6
    And the "fingerprint of Structure" reading of filter panel should be "Morgan"
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And no errors should have been logged
