@journey
Feature: A pane dragged out of the context panel becomes a viewer of the current molecule
  The panes of a molecule cell are computed by panel functions. One dragged out of the panel becomes an Info panel
  viewer of its own, floating over the table view and holding no table: it recalculates for each molecule that becomes
  current, and only for those, until it is closed. Clicking a cell of another semantic type or a column header leaves it
  on the last molecule. As a viewer it docks into the view and comes back with the layout, floating or docked. Out of
  sight it does not recalculate, and shows the latest molecule once it is back. Once closed, a later molecule does not
  bring it back.

  Here on the Properties pane, whose formula tells the molecules apart, in either Chemical Engine. The Gasteiger Partial
  Charges pane behaves the same; guides/following-panels.feature films it in Crux mode, as in RDKit mode it is a
  server-side Python script: see the bdd library's CLAUDE.md, "What never becomes a feature". The
  rules behind this (the function's condition and units, a result that comes late, Close All, undoing the close, a query
  or a converter behind the pane, a view that cannot hold viewers) are the platform's tests
  `xamgle | Context panel | Dragged out pane`.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles-50 dataset
    And the context panel is open
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    Then row 1 should be current
    And the context panel should show the current cell
    When user expands Chemistry accordion header in context panel
    And user drags Properties accordion header in context panel out of the accordion
    Then "Properties" floating panel should contain the text "C24H23NO3"

  Scenario: The dragged-out pane recalculates for each molecule
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And the context panel should show the current cell
    And "Properties" floating panel should contain the text "C12H14"
    And "Properties" floating panel should not contain the text "C24H23NO3"
    When user clicks on the "cell 6 of canonical_smiles" area of grid
    Then row 6 should be current
    And the context panel should show the current cell
    And "Properties" floating panel should contain the text "C27H35NO3"
    And no errors should have been logged

  Scenario: Anything else made current has no effect on it
    When user clicks on the "cell 6 of canonical_smiles" area of grid
    Then row 6 should be current
    And "Properties" floating panel should contain the text "C27H35NO3"
    When user clicks on the "cell 3 of molregno" area of grid
    Then row 3 should be current
    And the context panel should show the current cell
    And "Properties" floating panel should contain the text "C27H35NO3"
    When user clicks on the "header molregno" area of grid
    Then the context panel should show "molregno"
    And "Properties" floating panel should contain the text "C27H35NO3"
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    Then row 1 should be current
    And the context panel should show the current cell
    And "Properties" floating panel should contain the text "C24H23NO3"
    And no errors should have been logged

  Scenario: It comes back with the layout while it floats
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And "Properties" floating panel should contain the text "C12H14"
    When user saves the layout of the current table view
    And user loads the saved layout
    Then "Properties" floating panel should contain the text "C12H14"
    When user clicks on the "cell 6 of canonical_smiles" area of grid
    Then row 6 should be current
    And the context panel should show the current cell
    And "Properties" floating panel should contain the text "C27H35NO3"
    And no errors should have been logged

  Scenario: It docks into the view, follows there, and comes back with the layout docked
    When user docks info panel viewer to the right edge of the view
    Then "Properties" floating panel should be absent
    And "Properties" dock panel should contain the text "C27H35NO3"
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And the context panel should show the current cell
    And "Properties" dock panel should contain the text "C12H14"
    When user saves the layout of the current table view
    And user loads the saved layout
    Then "Properties" floating panel should be absent
    And "Properties" dock panel should contain the text "C12H14"
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    Then row 1 should be current
    And the context panel should show the current cell
    And "Properties" dock panel should contain the text "C24H23NO3"
    And no errors should have been logged

  Scenario: Out of sight it waits, and shows the latest molecule once it is back
    Given user opens smiles-50 dataset keeping the first 10 rows as "more molecules"
    And the context panel is open
    When user clicks on the "cell 6 of canonical_smiles" area of grid
    Then row 6 should be current
    And the context panel should show the current cell
    When user switches to the "smiles-50" table view
    Then "Properties" dock panel should contain the text "C27H35NO3"
    And no errors should have been logged

  Scenario: Closing it ends it, and nothing keeps it
    When user hovers over the text "Molecular formula" in "Properties" dock panel
    And user closes "Properties" dock panel, watching that nothing keeps it
    Then info panel viewer should be absent
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And the context panel should show the current cell
    And info panel viewer should be absent
    And nothing should keep the closed panel
    And no errors should have been logged
