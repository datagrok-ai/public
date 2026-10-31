@journey @realizes:filters.cp.chem-and-bio-filters @serial
Feature: The substructure card with the sketcher, other cards and a reopened panel
  On spgi-100 with an empty filter panel and a Series card, the Structure header dropped on the panel
  adds an empty substructure card on top and filters nothing. With Filter as you draw cleared, a
  molecule typed into the card's sketcher reaches the card only on OK: pyridine keeps 17 rows, and
  piperidine typed over it keeps them until OK brings 15. With the option checked, pyridine reaches
  the grid while the sketcher is still open. With pyridine applied, Align and Highlight cleared in the
  sketcher and the grid still, checking both repaints the grid while the sketcher is open, and the
  same 17 rows pass (TestTrack Chem/chem-filters.md, step 9, claim 1). Pyridine and the Pyrrolidines category together keep the
  5 rows that satisfy both; switching the substructure card off leaves the 21 Pyrrolidines rows and
  switching it on brings the 5 back. A switched-off card survives closing and reopening the panel:
  still off, still holding pyridine, and on again it keeps the 5 rows. Cancel in the sketcher with
  Filter as you draw on gives back what the card held: nothing for a new structure, pyridine for an
  edit (public#697). The panel's master switch off lets every row through and keeps the card's
  structure, and on again brings the 5 rows back (public#3328).
  Serial (`@serial`): Filter as you draw is the account's setting (`chem-filter` user settings), and this feature clears
  it; the features that filter while the sketcher is open take turns with it.

  Background:
    Given user is logged in
    And the "chem-filter" user settings are put back at feature end
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

  @realizes:public#697
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

  @realizes:public#697
  Scenario: Cancel on an edited structure drawn with Filter as you draw gives the old structure back (public#697)
    When user clicks on the "card Structure" area of filter panel
    Then sketcher dialog should be visible
    And "Filter as you draw" input in sketcher dialog should be checked
    When user clears molecule input of sketcher dialog
    And user types "C1CCNCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then 15 rows should pass the filter
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be hidden
    And 17 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And no errors should have been logged

  Scenario: Align and Highlight turned on with the sketcher open repaint the grid and keep the rows
    When user clicks on the "card Structure" area of filter panel
    Then sketcher dialog should be visible
    When user unchecks "Align" input in sketcher dialog
    And user unchecks "Highlight" input in sketcher dialog
    Then "Align" input in sketcher dialog should not be checked
    And "Highlight" input in sketcher dialog should not be checked
    And 17 rows should pass the filter
    When user takes a snapshot of grid
    Then grid should not have repainted
    When user checks "Align" input in sketcher dialog
    And user checks "Highlight" input in sketcher dialog
    Then grid should have repainted
    And sketcher dialog should be visible
    And 17 rows should pass the filter
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

  @realizes:public#3328
  Scenario: The panel switched off and on again filters by the structure the card holds (public#3328)
    When user hovers over filter panel
    And user unchecks master of filter panel
    Then the "active" reading of filter panel should be "false"
    And all rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    When user hovers over filter panel
    And user checks master of filter panel
    Then the "active" reading of filter panel should be "true"
    And 5 rows should pass the filter
    And the "filtering of Structure" reading of filter panel should be "true"
    And no errors should have been logged
