Feature: Rows and cells of a Hit Design campaign
  The campaign's grid is the platform's, with the package's row actions on top: a new row (the ribbon's +, or the
  Molecule cell's context menu) gets the template's first stage wherever the Stage is empty and opens the sketcher on
  its empty molecule; a Molecule cell's menu adds Add new row, Duplicate molecule and, for a cell with a molecule,
  Re-Run Calculations, which recomputes that row with the template's compute (here the JS-only Chemical Properties,
  MW) and saves. The Stage cell, a choice of the template's stages, is stage-choices.feature.

  Drawing the molecule of a new row, and the V-iD the drawing registers, is the UsageAnalysis feature
  apps/hit-design-crux.feature; it is not repeated here. Not claimed, as a known bug: Duplicate molecule gives the copy
  no V-iD and no computed values and pops the sketcher up over the filled row.

  The campaign is made for the feature through the JS API and opened by its link (bindings/campaign.ts).

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Grid" with the stages "Idea, Design, Make" and a table of:
      | Molecule | Stage  | V-iD | MW    |
      | CCO      | Design |      | 46.07 |
      | c1ccccc1 | Design |      |       |
    When user opens the campaign "BDD Grid" by its link
    Then the table should have 2 rows

  Scenario: Add new row gives the new row the first stage and opens the sketcher on its empty molecule
    When user clicks on "Add new row" icon in toolbar
    Then sketcher dialog should be visible
    And the table should have 3 rows
    And the value of "Stage" column in row 3 should be "Idea"
    And the value of "Stage" column in row 1 should be "Design"
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be absent
    And the table should have 3 rows

  Scenario: A Molecule cell's menu adds the package's row actions, and another column's menu does not
    When user right-clicks on the "cell 1 of Molecule" area of grid
    Then the open menu should list "Add new row"
    And the open menu should list "Duplicate molecule"
    And the open menu should list "Re-Run Calculations"
    When user closes the context menu
    And user right-clicks on the "cell 1 of Stage" area of grid
    Then the open menu should list "Sort..."
    And the open menu should not list "Duplicate molecule"
    And the open menu should not list "Re-Run Calculations"

  Scenario: Duplicate molecule appends the cell's molecule as a new row in the first stage
    When user picks "Duplicate molecule" from the context menu of the "cell 1 of Molecule" area of grid
    Then the table should have 3 rows
    And the value of "Molecule" column in row 3 should be "CCO"
    And the value of "Stage" column in row 3 should be "Idea"

  Scenario: Re-Run Calculations computes the template's property for that row alone and saves it
    Then the value of "MW" column in row 2 should be ""
    When user picks "Re-Run Calculations" from the context menu of the "cell 2 of Molecule" area of grid
    Then the campaign "BDD Grid" should be saved with a value matching "^78[.]1" in row 2 of the "MW" column
    And the value of "MW" column in row 2 should match "^78[.]1"
    And the value of "MW" column in row 1 should match "^46[.]0[67]"
    And no error or warning balloon should have been shown
