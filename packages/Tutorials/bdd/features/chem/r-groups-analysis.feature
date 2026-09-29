@tutorials @serial @realizes:tutorials.r-groups-analysis
Feature: The R-Groups Analysis tutorial
  Walks Cheminformatics > R-Groups Analysis from its card to the end: the analysis from the Chem menu
  with the scaffold found by MCS, the trellis of pie charts it adds set up from the context panel,
  a segment clicked, rows Shift-dragged in the grid and their distributions explored, the trellis
  switched to histograms and re-split. Each step is claimed as ticked and as done — the R-group
  columns and the trellis, the inner Category and Value, the selection, the Distributions pane,
  the inner viewer type, the X axis.
  Translated from playwright-tests/e2e/tutorials/r-groups-analysis.test.ts, which clicked across the
  trellis until a segment happened to be hit and forced the row group current through the status bar.
  MCS and the decomposition run in the browser (RDKit).

  Fixed in the tutorial for this translation: the counter never reached its last step (14 declared
  for 15 actions); the gear step completed on a click on any gear of the page; the Category and
  Value steps watched the context panel's DOM for a text and now read the trellis's inner viewer;
  the Pie chart tab's hint was captured when its step began; a stale Jupyter prerequisite (both
  computations run in the browser) is gone; the hover step bound its listener to a pane that did
  not exist yet when the step began, threw, and slept a second — it now follows the pointer and the
  tooltip (and "dictibutions" is spelt right).
  Needs the trellis plot's positional `cell body <column>,<row>` area (core): the categories here are
  SMILES, which no phrase can spell.
  The context panel drops a new current object within 2 s of a property edit (by design, GROK-21024),
  which a feature is faster than: a grid cell click, which releases that guard, comes before the
  Shift+Drag and before the gear that follows the Histogram pick.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "R-Groups Analysis" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the R-Groups Analysis tutorial
    When user starts the "R-Groups Analysis" tutorial
    Then the tutorial progress should be 1 of 15
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then the tutorial step "On the Top Menu, click Chem > Analyze > R-Groups Analysis..." should be done
    And "R-Groups Analysis" dialog should be visible
    When user clicks on "MCS" button in "R-Groups Analysis" dialog
    Then the tutorial step "Click MCS" should be done
    # the step ticks on the click; the scaffold lands in the sketcher when MCS ends, and OK before that
    # only says "No core was provided"
    And "R-Groups Analysis" dialog should have finished updating
    When user clicks on OK button in "R-Groups Analysis" dialog
    Then the tutorial step "Click OK" should be done
    And the tutorial step "Wait for the analysis to complete" should be done
    And trellis plot viewer should be visible
    And the table should have a column "R1"
    And the "inner viewer type" reading of trellis plot viewer should be "Pie chart"

    When user hovers over trellis plot viewer
    And user clicks on settings icon of trellis plot viewer
    Then the tutorial step "In the trellis plot, click the gear icon for the embedded viewer" should be done
    And the context panel should show "Trellis plot"
    When user clicks on "Pie chart" tab in context panel
    Then the tutorial step "Go to Pie chart tab" should be done
    When user selects "LC/MS" in "Category" property
    Then the tutorial step "Under Pie chart tab > Data, set Category to LC/MS" should be done
    And "categoryColumnName" inner property of trellis plot viewer should be "LC/MS"

    When user clicks on the "cell body 1,1" area of trellis plot viewer
    Then the tutorial step "Click any segment on a pie chart" should be done
    And some rows should be selected
    When user presses Escape
    Then the tutorial step "Press Escape" should be done

    When user clicks on the "cell 1 of R1" area of grid
    Then the context panel should show the current cell
    When user drags from the "row header 1" area to the "row header 7" area of grid holding Shift
    Then the tutorial step "In the grid, press Shift+Drag Mouse Down" should be done
    And 7 rows should be selected
    When user expands "Distributions" pane in context panel
    Then the tutorial step "On the Context Panel, expand the Distributions pane" should be done
    When user hovers over "Distributions" pane in context panel
    Then the tutorial step "In the pane, hover over line charts to see distributions" should be done
    And tooltip should be visible

    When user picks "Histogram" in the viewer selector of trellis plot viewer
    Then the tutorial step "In the top-left corner of the trellis plot, select Histogram." should be done
    And the "inner viewer type" reading of trellis plot viewer should be "Histogram"
    When user clicks on the "cell 8 of R1" area of grid
    Then the context panel should show the current cell
    When user hovers over trellis plot viewer
    And user clicks on settings icon of trellis plot viewer
    And user clicks on "Histogram" tab in context panel
    And user selects "In-vivo Activity" in "Value" property
    Then the tutorial step "Set Value to In-Vivo Activity" should be done
    And "valueColumnName" inner property of trellis plot viewer should be "In-vivo Activity"

    When user picks "R4" in the column selector at the "x selector 1" area of trellis plot viewer
    Then the tutorial step "Set the value for the X axis to R4" should be done
    And "X Column Names" property of trellis plot viewer should be "R4"

    And the "R-Groups Analysis" tutorial should be completed
    And the tutorial should have listed 15 steps
    And the tutorial progress should be 15 of 15
    And no hint should be shown
    And no errors should have been logged
