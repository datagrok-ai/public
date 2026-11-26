@guide @tutorial
Feature: How-to · Curves
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Curves package's Assay_Curves table (43 dose-response curves,
  a 4PL fit each). A click on a curve shows Context Panel > Fit: the statistics of that curve, each with a "+" whose
  tooltip says "Calculate <name> for the whole column"; the "+" of IC50 and of R² add the calculated columns "Fitted
  Curve 1 ic50" and "Fitted Curve 1 rSquared" right after the curves, one value per row. Sorted by R² ascending, the
  poorest fits come first. Off camera, before the first step: the table opened, the grid showing only the molecule,
  the curves and the curve's name (the file's own AUC, Hill, Min, Max and IC50 columns are hidden, so the IC50 shown
  is the one the clip adds), the toolbox hidden, the context panel open. Nothing is made on the server.

  Scenario: IC50 and other statistics as columns
    Turn a curve's IC50 and R² into columns with a value for every curve, then sort to find the poorest fits.

    Given user is logged in
    And simple mode is off
    And the "Curves" package is installed
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens assay-curves dataset
    And the toolbox pane is hidden
    And the grid shows only the columns "Molecule, Fitted Curve, Compound|Assay|Target"
    And the context panel is open
    # caption: Click a curve, then open Fit in the context panel
    When user clicks on the "cell 1 of Fitted Curve" area of grid beside its curve
    Then "Fit" pane in context panel should be visible
    When user expands Fit accordion header in context panel
    Then "IC50" input in context panel should be visible
    When the video says "Fit: this curve's statistics, from its fitted parameters"
    # caption: The + next to a value adds it for every curve
    And user hovers over first "Calculate ic50 for the whole column" icon in context panel
    And user clicks on first "Calculate ic50 for the whole column" icon in context panel
    Then the table should have a column "Fitted Curve 1 ic50"
    And "Fitted Curve 1 ic50" column should have type "double"
    And "Fitted Curve 1 ic50" column should have at least 40 distinct values
    When the video says "A new column with every curve's IC50, too small to read yet"
    # caption: R² tells how well each curve fits its points
    And user clicks on first "Calculate rSquared for the whole column" icon in context panel
    Then the table should have a column "Fitted Curve 1 rSquared"
    # caption: IC50s in molar are tiny: show them in scientific notation
    When user picks "Format > scientific (1.63E-7)" from the context menu of the "header Fitted Curve 1 ic50" area of grid
    And user closes the context menu
    Then the "Fitted Curve 1 ic50" cell of row 1 should be displayed as "1.63E-7"
    # caption: Sort by R²: the poorest fits come first
    When user picks "Sort > Ascending" from the context menu of the "header Fitted Curve 1 rSquared" area of grid
    And user closes the context menu
    Then the "sort direction" reading of grid should be "ascending"
    # caption: R² near 0: the curve does not fit these points
    # (cells are named by table row: after the sort, row 7, R² 0.00, is on top)
    When user hovers over the "cell 7 of Fitted Curve 1 rSquared" area of grid
    Then no errors should have been logged
