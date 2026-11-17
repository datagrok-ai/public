@tutorials @serial @realizes:apps.tutorials
Feature: The Tutorials application
  The app as a whole, as the TestTrack case "Tutorials" (Apps) walks it: launched from Browse > Apps its
  view lists every track and every tutorial; a tutorial runs in a bar on top, with all its steps in a
  panel on demand, and closing it returns to the view; a tutorial finished to its end offers the next
  one of its track, which starts at its first step; the finished one's card says it is done. In the
  Browse tree the app lists its tracks and tutorials, and selecting one previews it without starting it.
  Each tutorial walked to its end, with every step claimed, is a feature of its own (bdd/features/*) —
  the case's "Complete all tutorials", "no steps are skipped or missing" and "the highlighted elements"
  are claimed there, step by step and hint by hint.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "recentViewerSettings" user settings are put back at feature end

  Scenario: Browse > Apps > Tutorials lists every track and tutorial
    Given the Tutorials app is closed
    When user clicks on Apps tree node inside browse tree
    And user double-clicks on Tutorials gallery card
    Then Tutorials view should be visible
    And Tutorials panel should be absent
    And there should be 7 visible tutorial tracks
    And there should be 20 visible tutorial cards
    And "Exploratory data analysis" tutorial track should be visible
    And "Scatter Plot" tutorial card should be visible
    And no errors should have been logged

  Scenario: A running tutorial shows its current step in a bar
    Given the "Scatter Plot" tutorial is not completed yet
    And the Tutorials app is open
    When user clicks on "Scatter Plot" tutorial card
    Then tutorial bar should contain text "Open scatter plot"
    And Tutorials panel should be absent
    And the tutorial progress should be 1 of 10
    When user clicks on all steps icon
    Then Tutorials panel should be visible
    And tutorial bar should be above Tutorials panel
    When user clicks on all steps icon
    Then Tutorials panel should be absent
    And the tutorial step "Open scatter plot" should be current
    When user closes the tutorial
    Then tutorial bar should be absent
    And Tutorials panel should be absent
    And Tutorials view should be visible
    And no errors should have been logged

  Scenario: A finished tutorial offers the next one of its track, which starts
    Given the "Scatter Plot" tutorial is not completed yet
    And the "Embedded Viewers" tutorial is not completed yet
    And the Tutorials app is open
    When user starts the "Scatter Plot" tutorial
    And user clicks on scatter-plot icon in toolbox
    And user picks "HEIGHT" in the "x" column selector of scatter plot viewer
    And user picks "WEIGHT" in the "y" column selector of scatter plot viewer
    And user picks "AGE" in the "size" column selector of scatter plot viewer
    And user picks "SEX" in the "color" column selector of scatter plot viewer
    And user drags a zoom box over the "view" area of scatter plot viewer
    And user double-clicks on empty plot space of scatter plot viewer
    And user clicks on the "marker of row 11" area of scatter plot viewer
    And user drags a selection box over the "view" area of scatter plot viewer
    And user presses Escape in scatter plot viewer
    Then the "Scatter Plot" tutorial should be completed
    And tutorial bar should contain text "Next \"Embedded Viewers\""
    # the next is the first tutorial of the track not completed yet, after the finished one first
    When user clicks on Start button in tutorial bar
    Then tutorial title should contain text "Embedded Viewers"
    And the tutorial progress should be 1 of 9
    And the tutorial step "Open scatter plot" should not be done yet
    When user closes the tutorial
    Then the "Scatter Plot" tutorial card should show it is done
    And no errors should have been logged

  Scenario: Selecting a tutorial in the Browse tree previews it and starts nothing
    Given the Tutorials app is closed
    And Apps tree node inside browse tree is expanded
    And Apps---Tutorials tree node inside browse tree is expanded
    And Apps---Tutorials---Exploratory-Data-Analysis tree node inside browse tree is expanded
    When user clicks on Apps---Tutorials---Exploratory-Data-Analysis---Scatter-Plot tree node inside browse tree
    Then "Start" button should be visible
    And tutorial bar should be absent
    When user clicks on Apps---Tutorials---Exploratory-Data-Analysis tree node inside browse tree
    Then Tutorials view should be visible
    And there should be 1 visible tutorial tracks
    And tutorial bar should be absent
    And no errors should have been logged
