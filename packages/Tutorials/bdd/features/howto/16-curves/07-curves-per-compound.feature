@guide @tutorial
Feature: How-to · Curves
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Curves package's Assay_Curves table (43 curves of 26
  compounds, one to four runs each). A trellis plot from the toolbox, "Curves" (the Multi curve viewer) picked as the
  chart in every cell and the cells split by Compound ID across (the split down, Assay Name by default, removed with
  Reset Y columns in its selector's context menu), overlays each compound's runs in a cell of its own. In the
  settings, Misc > On Click = Filter makes a click on a cell filter the table to that compound: compound_id 20 has
  four runs (rows 34 to 37). Off camera, before the first step: the table opened, the toolbox shown, Browse and the
  context panel closed. Nothing is made on the server.

  Scenario: Curves per compound in a trellis plot
    One small chart per compound, each overlaying that compound's runs; click a chart to filter the table to it.

    Given user is logged in
    And simple mode is off
    And the "Curves" package is installed
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens assay-curves dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A trellis plot draws one small chart per group"
    And user clicks on trellis plot icon on toolbox
    Then trellis plot viewer should be added to the open tableview
    # caption: Hide the toolbox to give the plot room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Curves: the chart drawn in every cell
    When user picks "Curves" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "MultiCurveViewer"
    # caption: Split by compound across
    When user picks "Compound ID" in the column selector at the "x selector 1" area of trellis plot viewer
    # caption: Drop the split down: right-click, Reset Y columns
    And user picks "Reset Y columns" from the context menu of the "y selector 1" area of trellis plot viewer
    And user closes the context menu
    Then the "cells" reading of trellis plot viewer should be at least 20
    And trellis plot viewer should have a "cell compound_id 20" area
    When the video says "Each compound's runs, overlaid in a cell of its own"
    # caption: Settings: the gear on the plot's title bar
    When user hovers over trellis plot viewer
    And user clicks on settings icon of trellis plot viewer
    # caption: Misc > On Click: Filter
    And user expands "Misc" category in context panel
    And user selects "Filter" in "On Click" property in context panel
    # caption: Click a compound's chart: the table shows only its runs
    And user clicks on the "cell body compound_id 20" area of trellis plot viewer
    Then 4 rows should pass the filter
    And all rows where "Compound ID" is "compound_id 20" should pass the filter
    And no errors should have been logged
