@guide @tutorial
Feature: How-to · Curves
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Curves package's sample Assay_Curves.csv, opened from Browse
  (Files > App Data > Curves > demo), holds 43 dose-response curves, one per compound, assay, target and run: the
  "Fitted Curve" column is recognized as curves (semantic type fit) and each cell draws its points and a 4PL fit on a
  log concentration axis. The pointer on a point shows its concentration and response; a click on the cell (away
  from its points: a click on a point would leave it out of the fit) shows the curve in the context panel (Chart, and
  Fit with the statistics fitted in the browser); a double-click opens the chart in the Edit chart dialog. Off camera,
  before the first step: Browse opened down to App Data, the toolbox hidden, the context panel closed. Nothing is
  made on the server.

  Scenario: Dose-response curves in a table
    Open a table of fitted curves, read a single measurement, and look at one curve and its statistics up close.

    Given user is logged in
    And simple mode is off
    And the "Curves" package is installed
    And the "Chem" package is installed
    And the package autostarts have completed
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    # caption: Sample curves: Files > App Data > Curves > demo
    When user expands Files---App-Data---Curves tree node inside browse tree
    And user expands Files---App-Data---Curves---demo tree node inside browse tree
    And user double-clicks on Files---App-Data---Curves---demo---Assay_Curves.csv tree node inside browse tree
    Then the "Assay_Curves" view should be current
    And the table should have 43 rows
    And "Fitted Curve" column should have semantic type "fit"
    When the video says "Each cell draws a dose-response curve: the points and the fit"
    # caption: Point at a measurement: its concentration and response
    And user hovers over point 12 of the curve in the "cell 1 of Fitted Curve" area of grid
    # caption: Click a curve, then open the context panel
    And user clicks on the "cell 1 of Fitted Curve" area of grid
    Then no point of the curve in row 1 of the "Fitted Curve" column should be left out
    When user clicks on "Context Panel" status bar toggle
    Then "Chart" pane in context panel should be visible
    # caption: Chart: the same curve, larger, with its axes labelled
    When user expands Chart accordion header in context panel
    # caption: Fit: IC50, Hill slope, R² and more, fitted in your browser
    And user expands Fit accordion header in context panel
    Then "IC50" input in context panel should be visible
    # caption: Double-click a curve to open it in a window of its own
    When user double-clicks on the "cell 1 of Fitted Curve" area of grid
    Then "Edit chart" dialog should be visible
    And no point of the curve in row 1 of the "Fitted Curve" column should be left out
    When user clicks on "Close" button in "Edit chart" dialog
    Then "Edit chart" dialog should be absent
    And no errors should have been logged
