@guide @tutorial
Feature: How-to · Molecules
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). R-group analysis: Chem > Analyze > R-Groups Analysis...
  opens a dialog with a sketcher for the core; its MCS button fills in the largest substructure the molecules
  share, and OK splits every molecule into that core and its R-groups: one column per R-group and a trellis plot
  of the R-group combinations (RDKit's decomposition, in the browser). sar_small (Chem's demo of 200 molecules of
  one series) is opened off camera with the grid showing smiles, ID and In-vivo Activity; the new columns live in
  the page's table only and nothing is made on the server.

  Scenario: Analyze R-groups
    Find the common core of a series and split each molecule into the core and its R-groups.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens sar_small dataset
    And the grid shows only the columns "smiles, ID, In-vivo Activity"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Chem > Analyze > R-Groups Analysis
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then "R-Groups Analysis" dialog should be visible
    # caption: Draw the core, or let MCS find the common one
    When user clicks on MCS button in "R-Groups Analysis" dialog
    Then "R-Groups Analysis" dialog should have finished updating
    When the video says "MCS: the largest substructure all the molecules share" at the top
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then a new column matching "^R1" should have been added
    And trellis plot viewer should be visible
    When the video says "One column per R-group, and a plot of how they combine"
    Then the table should have 200 rows
    And no errors should have been logged
