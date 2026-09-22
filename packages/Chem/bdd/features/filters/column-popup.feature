@journey @realizes:filters.cp.chem-and-bio-filters
Feature: A column popup's substructure filter moved to the filter panel
  On spgi-100 with an empty filter panel, the Column options icon of the Structure header opens a
  popup titled Structure with a substructure filter of its own. Pyridine typed there keeps 17 rows
  while the panel still has no Structure card. Add filter moves the criterion to the panel: the
  popup's filter goes, the panel gets a Structure card holding pyridine, and the 17 rows stay. With
  the popup dismissed, switching the panel card off lets every row through, so no filter is left
  behind in the popup, and switching it on brings the 17 rows back (GROK-14952).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens spgi-100 dataset

  Scenario: The popup's filter narrows the rows on its own
    When user opens an empty filter panel
    And user hovers over the "header Structure" area of grid
    And user clicks on "Column options" icon in grid
    Then column popup should be visible
    And title of column popup should have text "Structure"
    And "Structure" filter card should be absent
    When user types "c1ccncc1" into molecule input of column popup
    And user presses Enter in molecule input of column popup
    Then 17 rows should pass the filter
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccncc1"
    And "Structure" filter card should be absent
    And no errors should have been logged

  Scenario: Add filter moves the criterion to the panel
    When user clicks on "Add filter" action in column popup
    Then "Structure" filter card should be visible
    And molecule input of column popup should be absent
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And 17 rows should pass the filter
    When user presses Escape
    Then column popup should be absent
    And 17 rows should pass the filter
    And no errors should have been logged

  Scenario: The panel card switched off lets every row through
    When user hovers over "Structure" filter card
    And user unchecks checkbox of "Structure" filter card
    Then "Structure" filter card should be disabled
    And all rows should pass the filter
    When user hovers over "Structure" filter card
    And user clicks on checkbox of "Structure" filter card
    Then 17 rows should pass the filter
    And no errors should have been logged
