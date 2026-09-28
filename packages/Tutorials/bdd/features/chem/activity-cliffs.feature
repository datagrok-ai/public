@tutorials @serial @realizes:tutorials.activity-cliffs
Feature: The Activity Cliffs tutorial
  Walks Cheminformatics > Activity Cliffs from its card to the end: the analysis from the Chem menu
  with its default settings, the plot it draws explored — a molecule hovered, the plot narrowed to
  the cliffs, zoomed, a cliff line hovered and clicked — the pair opened from the context panel, and
  the table of cliffs docked, enlarged and clicked. Each step is claimed as ticked and as done — the
  plot and its 15 cliffs, the tooltip, the filter, the zoom, the line under the pointer, the pair's
  pane, the current row, the table and its height.
  Translated from playwright-tests/e2e/tutorials/activity-cliffs.test.ts, which looked for the green
  lines by scanning the canvas for their colour and swept the plot until a tooltip came up. The
  lines renderer reports each line it drew as an area; the analysis runs in the browser.

  Fixed in the tutorial for this translation: the Show only cliffs step completed on a click on any
  switch of the page and its hint was captured when the step began; it now waits for the plot's
  own filter.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Activity Cliffs" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Activity Cliffs tutorial
    When user starts the "Activity Cliffs" tutorial
    Then the tutorial progress should be 1 of 12
    When user picks "Chem > Analyze > Activity Cliffs..." from the top menu
    Then the tutorial step "On the Top Menu, click Chem > Analyze > Activity Cliffs..." should be done
    And "Activity Cliffs" dialog should be visible
    When user clicks on OK button in "Activity Cliffs" dialog
    Then the tutorial step "Click OK" should be done
    And the tutorial step "Wait for analysis to complete" should be done
    And scatter plot viewer should be visible
    And the "cliffs" reading of scatter plot viewer should be 15

    When user hovers over the first "marker" area of scatter plot viewer
    Then the tutorial step "Hover over data points for molecule information" should be done
    And tooltip should be visible
    When user toggles "Show only cliffs" input in scatter plot viewer
    Then the tutorial step "To view only the cliffs, toggle Show only cliffs." should be done
    And the "only cliffs" reading of scatter plot viewer should be "true"

    When user drags a zoom box over the "view" area of scatter plot viewer
    Then the tutorial step "Press Use Alt + Mouse Drag to zoom in" should be done

    When user hovers over the first free "line" area of scatter plot viewer
    Then the tutorial step "Hover over the green line to see the pair of molecules" should be done
    And tooltip should be visible
    When user clicks on that area of scatter plot viewer
    Then the tutorial step "Click on the green line connecting that molecule pair" should be done
    And "Cliff Details" pane in context panel should be visible
    # the line click made the pair's first molecule the current row already, and the step waits for the
    # current row to change: the other molecule is the one that moves it
    When user clicks on second cliff molecule in context panel
    Then the tutorial step "On the Context Panel, click any molecule" should be done

    When user clicks on "15 cliffs" button in scatter plot viewer
    Then the tutorial step "At the top right corner of the scatterplot, click 15 CLIFFS" should be done
    And the open tableview should have 2 grid viewers
    When user drags the top border of second grid viewer by 150 pixels up
    Then the tutorial step "Drag the top border of that table upwards to create more space for it" should be done
    When user clicks on the first "cell" area of second grid viewer
    Then the tutorial step "In the cliffs table, click any cell in the first row" should be done

    And the "Activity Cliffs" tutorial should be completed
    And the tutorial should have listed 12 steps
    And the tutorial progress should be 12 of 12
    And no hint should be shown
    And no errors should have been logged
