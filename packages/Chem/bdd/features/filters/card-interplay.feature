@journey @realizes:filters.cp.chem-and-bio-filters
Feature: The substructure card with the sketcher, other cards and a reopened panel
  On spgi-100 with an empty filter panel and a Series card, the Structure header dropped on the panel
  adds an empty substructure card on top and filters nothing. With Filter as you draw cleared, a
  molecule typed into the card's sketcher reaches the card only on OK: pyridine keeps 17 rows, and
  piperidine typed over it keeps them until OK brings 15. With the option checked, pyridine reaches
  the grid while the sketcher is still open. Pyridine and the Pyrrolidines category together keep the
  5 rows that satisfy both; switching the substructure card off leaves the 21 Pyrrolidines rows and
  switching it on brings the 5 back. A switched-off card survives closing and reopening the panel:
  still off, still holding pyridine, and on again it keeps the 5 rows.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens spgi-100 dataset

  Scenario: The Structure header dropped on the panel adds an empty card on top
    When user opens an empty filter panel
    And user adds a card for "Series" to the filter panel
    And user drags the "header Structure" area of grid onto the "view" area of filter panel
    Then "Structure" filter card should be visible
    And the "cards" reading of filter panel should be "Structure, Series"
    And the "type of Structure" reading of filter panel should be "Chem:substructureFilter"
    And the "structure of Structure" reading of filter panel should be ""
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: With Filter as you draw cleared the sketch reaches the card on OK
    When user clicks on "Sketch" text in "Structure" filter card
    And user unchecks "Filter as you draw" input in sketcher dialog
    Then "Filter as you draw" input in sketcher dialog should not be checked
    When user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then sketcher dialog should be visible
    And the "structure of Structure" reading of filter panel should be ""
    And all rows should pass the filter
    When user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And no errors should have been logged

  Scenario: With Filter as you draw cleared a further edit waits for OK
    When user clicks on the "card Structure" area of filter panel
    Then sketcher dialog should be visible
    And "Filter as you draw" input in sketcher dialog should not be checked
    When user clears molecule input of sketcher dialog
    And user types "C1CCNCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And 17 rows should pass the filter
    When user clicks on OK button in sketcher dialog
    Then 15 rows should pass the filter
    And the filter should pass exactly the molecules of "Structure" column containing "C1CCNCC1"
    And no errors should have been logged

  Scenario: With Filter as you draw checked the edit reaches the grid with the sketcher open
    When user clicks on the "card Structure" area of filter panel
    And user checks "Filter as you draw" input in sketcher dialog
    Then "Filter as you draw" input in sketcher dialog should be checked
    When user clears molecule input of sketcher dialog
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then 17 rows should pass the filter
    And sketcher dialog should be visible
    When user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    And no errors should have been logged

  Scenario: The substructure and category cards intersect and switch off on their own
    When user clicks on the "category Pyrrolidines of Series" area of filter panel
    Then 5 rows should pass the filter
    And no rows where "Series" is "Triazoles" should pass the filter
    When user hovers over "Structure" filter card
    And user unchecks checkbox of "Structure" filter card
    Then "Structure" filter card should be disabled
    And 21 rows should pass the filter
    When user hovers over "Structure" filter card
    And user clicks on checkbox of "Structure" filter card
    Then 5 rows should pass the filter
    And no errors should have been logged

  Scenario: A switched-off card survives closing and reopening the panel
    When user hovers over "Structure" filter card
    And user unchecks checkbox of "Structure" filter card
    Then 21 rows should pass the filter
    When user clicks on close icon of filters viewer
    Then filter panel should be hidden
    When user clicks on filter icon in toolbar
    Then "Structure" filter card should be disabled
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And 21 rows should pass the filter
    When user hovers over "Structure" filter card
    And user clicks on checkbox of "Structure" filter card
    Then 5 rows should pass the filter
    And no errors should have been logged
