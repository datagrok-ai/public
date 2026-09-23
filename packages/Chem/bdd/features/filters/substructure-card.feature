@journey @realizes:filters.cp.chem-and-bio-filters
Feature: The substructure filter card and its search types
  On spgi-100 (100 rows, molecule column Structure) the filter panel opened from the toolbar has a
  substructure card on the molecule columns. Removing them with their close icons leaves no
  substructure card; Add Filter | Substructure Filter... with every column checked in its picker
  puts them back. Pyridine drawn on the Structure card keeps 17 rows and the header counter reads 1.
  The card offers seven search types: Contains and Not contains split the table into 17 and 83 rows,
  Included in and Not included in split it into 0 and 100, Exact, Stereo agnostic and Similar keep
  none (pyridine is no molecule of the table and is similar to none), and Contains again keeps the
  same 17 rows.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens spgi-100 dataset

  Scenario: The panel opens with substructure cards, and their close icons remove them
    When user clicks on filter icon in toolbar
    Then "Structure" filter card should be visible
    And the "type of Structure" reading of filter panel should be "Chem:substructureFilter"
    And the "type of Core" reading of filter panel should be "Chem:substructureFilter"
    When user hovers over "Structure" filter card
    And user clicks on close of "Structure" filter card
    And user hovers over "Core" filter card
    And user clicks on close of "Core" filter card
    And user hovers over "R1" filter card
    And user clicks on close of "R1" filter card
    And user hovers over "R2" filter card
    And user clicks on close of "R2" filter card
    And user hovers over "R3" filter card
    And user clicks on close of "R3" filter card
    And user hovers over "R100" filter card
    And user clicks on close of "R100" filter card
    And user hovers over "R101" filter card
    And user clicks on close of "R101" filter card
    Then the filter panel should have no "Chem:substructureFilter" filter card
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Add Filter | Substructure Filter... puts a card on every checked molecule column
    When user picks "Add Filter | Substructure Filter..." from the viewer menu of filter panel
    Then "Select columns..." dialog should be visible
    When user clicks on All label in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Structure" filter card should be visible
    And the "type of Structure" reading of filter panel should be "Chem:substructureFilter"
    And the "type of R101" reading of filter panel should be "Chem:substructureFilter"
    And the "structure of Structure" reading of filter panel should be ""
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Pyridine on the Structure card keeps the rows that contain it
    When user clicks on "Sketch" text in "Structure" filter card
    Then sketcher dialog should be visible
    When user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccncc1"
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: The search types split and nest the rows
    When user opens the settings of the "Structure" filter card
    Then the "Structure" filter card should offer search types "Contains, Included in, Exact, Stereo agnostic, Similar, Not contains, Not included in"
    When user picks search type "Not contains" in the "Structure" filter card
    Then 83 rows should pass the filter
    And the "search type of Structure" reading of filter panel should be "Not contains"
    When user picks search type "Included in" in the "Structure" filter card
    Then 0 rows should pass the filter
    When user picks search type "Not included in" in the "Structure" filter card
    Then all rows should pass the filter
    When user picks search type "Exact" in the "Structure" filter card
    Then 0 rows should pass the filter
    When user picks search type "Similar" in the "Structure" filter card
    Then 0 rows should pass the filter
    When user picks search type "Stereo agnostic" in the "Structure" filter card
    Then 0 rows should pass the filter
    When user picks search type "Contains" in the "Structure" filter card
    Then 17 rows should pass the filter
    And the "search type of Structure" reading of filter panel should be "Contains"
    And no errors should have been logged
