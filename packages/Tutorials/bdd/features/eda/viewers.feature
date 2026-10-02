@tutorials @serial @realizes:tutorials.viewers
Feature: The Viewers tutorial
  Walks Exploratory Data Analysis > Viewers from its card to the congratulations: the Add viewer
  gallery by tag and by search, four viewers from the toolbox, hover, selection and the current
  record shared between them, a scatter plot's properties, and its style cloned, picked up and
  applied to another plot. Each step is claimed as ticked and as done — the gallery filtered, the
  viewers added, the rows hovered, selected and made current, the property set on the first plot
  arriving on the third through Pick Up / Apply (not read back where it was written).
  Translated from playwright-tests/e2e/tutorials/viewers.test.ts, whose ribbon retries and event
  dispatched on a detached node were workarounds for tutorial bugs fixed since (GROK-20408).

  Fixed in the tutorial for this translation: three selection steps completed on any selection
  event, so the tail of one gesture could tick the next step — each waits now for a settled selection
  other than the one it began with, and the claims below hold the next step open meanwhile; the
  gallery's search box was captured when its step began.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "recentViewerSettings" user settings are put back at feature end
    And the package autostarts have completed
    And the "Viewers" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Viewers tutorial
    When user starts the "Viewers" tutorial
    Then the tutorial progress should be 1 of 20
    And "Add viewer" icon should be hinted
    When user clicks on "Add viewer" icon
    Then viewer gallery should be visible
    And the tutorial step "Click the Add viewer icon to open the gallery" should be done

    Then "Charts" viewer tag should be hinted
    When user remembers the number of visible viewer cards
    And user clicks on "Charts" viewer tag
    Then the tutorial step "Click the \"Charts\" tag to filter the viewers" should be done
    And there should be fewer visible viewer cards than remembered
    When user clicks on "Radar" viewer card
    Then the tutorial step "Select the Radar viewer" should be done
    And the current view should hold at least 2 viewers

    When user clicks on "Add viewer" icon
    Then the tutorial step "Open the viewer gallery again" should be done
    When user types "Sunburst" into viewer gallery search
    Then the tutorial step "Type \"Sunburst\" in the search box" should be done
    When user clicks on "Sunburst" viewer card
    Then the tutorial step "Select the Sunburst viewer" should be done
    And the open tableview should have 1 sunburst viewer

    When user clicks on scatter-plot icon in toolbox
    Then the tutorial step "Open scatter plot" should be done
    When user clicks on histogram icon in toolbox
    Then the tutorial step "Open histogram" should be done
    When user clicks on pie-chart icon in toolbox
    Then the tutorial step "Open pie chart" should be done

    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tutorial step "Hover over the histogram bins or scatter plot points" should be done
    And the "hovered row" reading of scatter plot viewer should not be "0"

    When user drags a selection box over the "view" area of scatter plot viewer
    Then the tutorial step "Select points on the scatter plot" should be done
    And some rows should be selected
    # the drag's last selection events must not tick the next step
    And the tutorial step "Select one of the bins on the histogram" should not be done yet
    When user remembers the "rows selected" reading of scatter plot viewer
    And user clicks on the "bin 3" area of histogram viewer
    Then the tutorial step "Select one of the bins on the histogram" should be done
    And the "rows selected" reading of scatter plot viewer should not be as remembered
    And the tutorial step "Click a Sunburst segment to select its rows" should not be done yet
    When user remembers the "rows selected" reading of scatter plot viewer
    # the sunburst's rings are SEX, CONTROL and RACE; its inner "M" segment selects the men
    And user clicks on the "segment M" area of sunburst viewer
    Then the tutorial step "Click a Sunburst segment to select its rows" should be done
    And the "rows selected" reading of scatter plot viewer should not be as remembered
    And only rows where "SEX" is "M" should be selected

    When user clicks on the "marker of row 11" area of scatter plot viewer
    Then the tutorial step "Click on a point to set the current record" should be done
    And the "hovered row" reading of scatter plot viewer should be the current row

    When user picks "Properties..." from the viewer menu of scatter plot viewer
    Then the tutorial step "Open the scatter plot's properties" should be done
    And the context panel should show "Scatter plot"
    When user sets "markerDefaultSize" property of scatter plot viewer to "12"
    Then the tutorial step "Change a few visual properties, e.g., the background color or marker size" should be done

    When user picks "General > Clone" from the context menu of the "view" area of scatter plot viewer
    Then the tutorial step "Clone the scatter plot" should be done
    And the open tableview should have 2 scatter plot viewers
    When user picks "Pick Up / Apply > Pick Up" from the context menu of the "view" area of first scatter plot viewer
    Then the tutorial step "Pick up the scatter plot's style" should be done
    When user clicks on scatter-plot icon in toolbox
    Then the tutorial step "Open scatter plot" should be done 2 times
    And the open tableview should have 3 scatter plot viewers
    When user picks "Pick Up / Apply > Apply" from the context menu of the "view" area of third scatter plot viewer
    Then the tutorial step "Apply the style to the new viewer" should be done
    And "markerDefaultSize" property of third scatter plot viewer should be "12"

    And the "Viewers" tutorial should be completed
    And the tutorial should have listed 20 steps
    And the tutorial progress should be 20 of 20
    And no hint should be shown
    And no errors should have been logged
