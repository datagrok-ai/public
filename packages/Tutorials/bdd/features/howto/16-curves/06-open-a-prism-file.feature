@guide @tutorial
Feature: How-to · Curves
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Curves package reads GraphPad Prism files (.pzfx): its sample
  dose-response-4pl.pzfx (Files > App Data > Curves in Browse) holds one XY table, "Dose Response": ten
  concentrations, replicates for Compound A, B, C (partial) and D (with outliers), two of D's values marked Excluded
  in Prism. Opened with a double-click, the file shows as a view of its name: a row per XY table, its data sets fitted
  as curves in one "Fitted Curve" cell, the excluded values as outliers (crosses). A double-click on the cell (beside
  its curves) opens the Edit chart dialog, large enough for the legend; a click on it and the Context Panel toggle
  show each compound's statistics in Fit. Off camera, before the first step: Browse opened down to App Data, the
  toolbox hidden, the context panel closed. Nothing is made on the server.

  Scenario: Open a GraphPad Prism file
    Open a Prism .pzfx file: its XY tables come in as fitted curves, with statistics for every data set.

    Given user is logged in
    And simple mode is off
    And the "Curves" package is installed
    And the package autostarts have completed
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    # caption: Prism files (.pzfx) open like any other file
    When user expands Files---App-Data---Curves tree node inside browse tree
    And user double-clicks on Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree
    Then the "dose-response-4pl.pzfx" view should be current
    When the video says "Each XY table is a row, its data sets fitted as curves"
    # caption: Double-click the cell to see the curves larger
    And user double-clicks on the "cell 1 of Fitted Curve" area of grid beside its curve
    Then "Edit chart" dialog should be visible
    When the video says "One curve per compound; values excluded in Prism are crosses"
    And user clicks on "Close" button in "Edit chart" dialog
    Then "Edit chart" dialog should be absent
    # caption: Click the cell, then open the context panel
    When user clicks on the "cell 1 of Fitted Curve" area of grid beside its curve
    And user clicks on "Context Panel" status bar toggle
    Then "Fit" pane in context panel should be visible
    # caption: Fit: R², IC50, slope and more for every compound
    When user expands Fit accordion header in context panel
    Then "IC50" input in context panel should be visible
    And no errors should have been logged
