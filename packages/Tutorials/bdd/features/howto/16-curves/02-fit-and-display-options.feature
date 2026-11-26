@guide @tutorial
Feature: How-to · Curves
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Curves package's Assay_Curves table (43 dose-response curves,
  a 4PL fit each): a click on a curve cell shows Context Panel > Options, whose Level (Column by default) says the
  options below it are set for the whole column. Fit Function switches every curve to another model (linear) and
  back; Confidence Interval draws a band round each fit; Droplines add IC80 to IC50; Statistics prints the IC50 and R²
  on the plot. A click on another curve of the column shows the same options. Off camera, before the first step: the
  table opened, its curve cells made larger (the dropline labels and the statistics need room), the toolbox hidden,
  the context panel open. The options live in the open table's column tags only: nothing is made on the server.

  Scenario: Change how curves are fitted and drawn
    Pick the model fitted to the curves, then add a confidence band, IC80 drop lines and statistics to every plot.

    Given user is logged in
    And simple mode is off
    And the "Curves" package is installed
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens assay-curves dataset
    And the toolbox pane is hidden
    And the rows of the grid are 230 px tall
    And the "Fitted Curve" column of the grid is 340 px wide
    And the context panel is open
    # caption: Click a curve: its options are in the context panel
    When user clicks on the "cell 1 of Fitted Curve" area of grid
    Then no point of the curve in row 1 of the "Fitted Curve" column should be left out
    And "Options" pane in context panel should be visible
    When user expands Options accordion header in context panel
    Then "Level" input in context panel should have value "Column"
    When the video says "Level: Column, so each change reaches every curve in the column"
    # caption: Fit Function: the model fitted to the points
    And user selects "linear" in "Fit Function" input in context panel
    Then "Fit Function" input in context panel should have value "linear"
    # caption: Back to the four-parameter dose-response fit
    When user selects "4pl-dose-response" in "Fit Function" input in context panel
    # caption: Confidence Interval: a band around each fit
    And user checks "Confidence Interval" input in context panel
    # caption: Droplines: IC80 next to IC50
    And user ticks "IC80" in "Droplines" input in context panel
    # caption: Statistics: values printed on every plot
    And user ticks "ic50" in "Statistics" input in context panel
    And user ticks "rSquared" in "Statistics" input in context panel
    # caption: Another curve of the column has the same options
    And user clicks on the "cell 2 of Fitted Curve" area of grid beside its curve
    Then no point of the curve in row 2 of the "Fitted Curve" column should be left out
    And "Fit Function" input in context panel should have value "4pl-dose-response"
    And "Confidence Interval" input in context panel should be checked
    And "ic50" should be ticked in "Statistics" input in context panel
    And no errors should have been logged
