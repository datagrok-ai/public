@guide @tutorial
Feature: How-to · Curves
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Curves package's Assay_Curves table, row 24 (compound_id 15,
  Collapse Assay, run_id 5): ten points falling from 105-123 % to -5 %, but point 5 (1.11E-9, 83.2 %) sits well
  below its neighbours (123.3 % and 121.1 %). The curves allow Click To Toggle, so a click on the point marks it an
  outlier (stored in the cell as `outlier: true`), draws it as a red cross and refits the curve without it; the R²
  and IC50 columns, calculated from the curves, follow: R² 0.93 -> 0.98, IC50 4.96E-8 -> 4.25E-8 (read on the stand).
  A second click puts the point back. Off camera, before the first step: the table opened, the R² and IC50 of the
  curves added as columns the way the Fit pane's "+" adds them (clip 03 shows it), IC50 in scientific notation, the
  file's own statistic columns hidden, the curve cells made larger, the grid scrolled to row 24, the toolbox hidden,
  the context panel closed. The change lives in the open table only: nothing is made on the server.

  Scenario: Leave out an outlier
    Click a point that does not belong to leave it out of the fit; the curve, its R² and its IC50 follow.

    Given user is logged in
    And simple mode is off
    And the "Curves" package is installed
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens assay-curves dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the "ic50" statistic of the "Fitted Curve" column is a column of the table
    And the "rSquared" statistic of the "Fitted Curve" column is a column of the table
    And user sets the format of "Fitted Curve 1 ic50" column in the grid to "scientific"
    And the grid shows only the columns "Molecule, Fitted Curve, Fitted Curve 1 rSquared, Fitted Curve 1 ic50, Compound|Assay|Target"
    And the rows of the grid are 230 px tall
    And the "Fitted Curve" column of the grid is 340 px wide
    And the grid of the current view starts at row 24
    Then the "Fitted Curve 1 rSquared" cell of row 24 should be displayed as "0.93"
    And the "Fitted Curve 1 ic50" cell of row 24 should be displayed as "4.96E-8"
    When the video says "R² and IC50 of every curve, calculated from its fit"
    # caption: R² 0.93: this curve misses one of its points
    And user hovers over the "cell 24 of Fitted Curve 1 rSquared" area of grid
    # caption: This point sits well below its neighbours
    And user hovers over point 5 of the curve in the "cell 24 of Fitted Curve" area of grid
    # caption: Click it to leave it out of the fit
    And user clicks on point 5 of the curve in the "cell 24 of Fitted Curve" area of grid
    Then point 5 of the curve in row 24 of the "Fitted Curve" column should be left out
    And the "Fitted Curve 1 rSquared" cell of row 24 should be displayed as "0.98"
    And the "Fitted Curve 1 ic50" cell of row 24 should be displayed as "4.25E-8"
    When the video says "Left out: a red cross, and the curve is fitted without it"
    # caption: R² follows: 0.98 now
    And user hovers over the "cell 24 of Fitted Curve 1 rSquared" area of grid
    # caption: So does IC50: from 4.96E-8 to 4.25E-8
    And user hovers over the "cell 24 of Fitted Curve 1 ic50" area of grid
    # caption: Click the cross to put the point back
    And user clicks on point 5 of the curve in the "cell 24 of Fitted Curve" area of grid
    Then point 5 of the curve in row 24 of the "Fitted Curve" column should be in the fit
    And the "Fitted Curve 1 rSquared" cell of row 24 should be displayed as "0.93"
    And the "Fitted Curve 1 ic50" cell of row 24 should be displayed as "4.96E-8"
    When the video says "Back in the fit: R² 0.93 and IC50 4.96E-8 again"
    Then no errors should have been logged
