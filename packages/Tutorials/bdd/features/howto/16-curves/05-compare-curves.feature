@guide @tutorial
Feature: How-to · Curves
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Curves package's Assay_Curves table, rows 34 to 37: the four
  runs (run_id 5, 4, 3, 2) of compound_id 20 in the Collapse Assay. Selected by their row headers, they are overlaid
  by the Multi curve viewer (Add viewer, search "curve"), which draws the curves of the selected rows, the current
  row and the row under the pointer, its legend naming each by its series (the run). Run 3 sits left of the other
  three (the file's own IC50 column: -8.54 against -7.98, -8.02 and -7.98). The pointer on row 38 (compound_id 21,
  run_id 5, drawn red) adds its curve; the clip ends with the pointer there (the curve stays after the pointer
  leaves the grid). Known wart, reported: the chart's title is the first overlaid cell's own title ("... | run_id 2")
  and the viewer's Title property is hidden in its settings, so it cannot be renamed.
  Off camera, before the first step: the table opened, the grid showing the molecule, the curves, Compound ID and Run
  ID, scrolled to row 34, the toolbox hidden, the context panel closed. Nothing is made on the server.

  Scenario: Compare curves in one chart
    Select the curves to compare and overlay them in one chart with the Multi curve viewer.

    Given user is logged in
    And simple mode is off
    And the "Curves" package is installed
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens assay-curves dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "Molecule, Fitted Curve, Compound ID, Run ID"
    And the grid of the current view starts at row 34
    # caption: Select the curves to compare: four runs of compound_id 20
    When user selects rows 34 to 37 by their row headers
    Then 4 rows should be selected
    And only rows where "Compound ID" is "compound_id 20" should be selected
    # caption: Add viewer, then search for curve
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "curve" into viewer gallery search
    # caption: Multi curve viewer: many curves in one chart
    And user clicks on first "MultiCurveViewer" card in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And the open tableview should have 1 MultiCurveViewer viewer
    And the "curves" reading of MultiCurveViewer viewer should be 4
    When the video says "The four runs overlaid, the legend names each run" at the bottom
    And the video says "Run 3 (squares) sits left of the others: a lower IC50" at the bottom
    # caption: Point at another row: its curve joins the comparison
    And user hovers over the "cell 38 of Compound ID" area of grid
    Then the "curves" reading of MultiCurveViewer viewer should be 5
    When the video says "In red: compound_id 21, against compound_id 20's runs" at the bottom
    Then no errors should have been logged
