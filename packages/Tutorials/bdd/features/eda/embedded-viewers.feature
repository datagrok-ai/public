@tutorials @serial @realizes:tutorials.embedded-viewers
Feature: The Embedded Viewers tutorial
  Walks Exploratory Data Analysis > Embedded Viewers from its card to the congratulations: three
  viewers from the toolbox, a scatter plot as the group tooltip of the others, a Trellis plot made
  from the pie chart with a scatter plot inside, and the inner plot's Color. Each step is claimed as
  ticked and as done on the platform — the tooltip really hosts the scatter plot and stops hosting
  it, the trellis really draws scatter plots, and the inner Color really repaints its cells.
  Translated from playwright-tests/e2e/tutorials/embedded-viewers.test.ts; the old spec checked the
  table tag, never the tooltip.

  Step 9 used to point at a gear icon no element carries (GROK-20419: no highlight, text about a
  control that is not there); the tutorial now points at the Color selector in the Trellis strip and
  says so, and the hint claims below hold it to that.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "recentViewerSettings" user settings are put back at feature end
    And the "Embedded Viewers" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Embedded Viewers tutorial
    When user starts the "Embedded Viewers" tutorial
    Then the tutorial progress should be 1 of 9
    When user clicks on scatter-plot icon in toolbox
    Then the tutorial step "Open scatter plot" should be done
    When user clicks on histogram icon in toolbox
    Then the tutorial step "Open histogram" should be done
    When user clicks on pie-chart icon in toolbox
    Then the tutorial step "Open pie chart" should be done
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 1 histogram viewer
    And the open tableview should have 1 pie chart viewer

    When user picks "Tooltip > Use as Group Tooltip" from the context menu of the "view" area of scatter plot viewer
    Then the tutorial step "Set the scatter plot as a tooltip viewer" should be done
    When user hovers over the "bin 3" area of histogram viewer
    Then the tutorial step "Hover over the histogram bins or pie chart segments" should be done
    And scatter plot viewer in tooltip should be visible

    When user picks "Tooltip > Remove Group Tooltip" from the context menu of the "view" area of scatter plot viewer
    Then the tutorial step "Reset the tooltip" should be done
    When user hovers over the "bin 4" area of histogram viewer
    Then tooltip should be visible
    And scatter plot viewer in tooltip should be absent

    When user picks "General > Use in Trellis" from the context menu of pie chart viewer
    Then the tutorial step "Open a Trellis plot from the pie chart's context menu" should be done
    And the open tableview should have 1 trellis plot viewer
    And the "inner viewer type" reading of trellis plot viewer should be "Pie chart"

    Then viewer selector in trellis plot viewer should be hinted
    When user picks "Scatter plot" in the viewer selector of trellis plot viewer
    Then the tutorial step "Set a scatter plot as an inner viewer" should be done
    And the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"

    # beside the Tutorials panel, the toolbox and the context panel the Trellis plot's dock is too
    # narrow for its strip: the inner Color selector lies past the dock's edge, hidden, and so is its
    # hint. The learner makes room first.
    Given the toolbox pane is hidden
    # F4 toggles the panel, which an earlier feature on the page may have left either way
    And the context panel is open
    When user presses F4
    Then context panel should be hidden
    And color column selector in trellis plot viewer should be hinted
    # the trellis splits by the pie chart's DIS_POP and by SEVERITY; its cells are viewers of their own,
    # so the inner Color is claimed on a cell's picture, which a numeric colour scale must change
    When user remembers the "cell signature RA | Critical" reading of trellis plot viewer
    And user picks "AGE" in the "color" column selector of trellis plot viewer
    Then the tutorial step "Set Color of the inner scatter plot to AGE" should be done
    And "colorColumnName" inner property of trellis plot viewer should be "AGE"
    And the "cell signature RA | Critical" reading of trellis plot viewer should not be as remembered

    And the "Embedded Viewers" tutorial should be completed
    And the tutorial should have listed 9 steps
    And the tutorial progress should be 9 of 9
    And no hint should be shown
    And no errors should have been logged
