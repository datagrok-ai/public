@guide @help:datagrok/solutions/domains/chem
Feature: Break a compound series into its core and R-groups
  A guide: the answer to "how do I split my series into the common core and its substituents, and
  see which combinations of substituents I have?". Chem > Analyze > R-Groups Analysis... opens a
  sketcher for the core, and its MCS button draws the maximum common substructure of the molecule
  column there, so the core need not be drawn by hand. OK decomposes every molecule on it: one new
  column per R-group position (R1, R2, ...), the R-groups highlighted in color in the structures,
  which are aligned on the core, and a trellis plot of R1 against R2 with a chart per combination.
  Demo: sar_small, 200 molecules of one series.

  Scenario: Find the core with MCS, then decompose the series on it
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens sar-small dataset
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then "R-Groups Analysis" dialog should be visible
    When user clicks on MCS button in "R-Groups Analysis" dialog
    Then "R-Groups Analysis" dialog should have finished updating
    When user clicks on OK button in "R-Groups Analysis" dialog
    Then a new column "R1" should have been added
    And a new column "R2" should have been added
    And trellis plot viewer should be visible
