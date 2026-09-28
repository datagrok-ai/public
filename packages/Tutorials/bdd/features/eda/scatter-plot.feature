@tutorials @serial @realizes:tutorials.scatter-plot
Feature: The Scatter Plot tutorial
  Walks the Exploratory Data Analysis > Scatter Plot tutorial the way a learner does, from its card
  to the congratulations. Every step is claimed twice: the tutorial ticks it (its entry checked, the
  progress moving) and the platform really did what the step asked — the column bound, the viewport
  zoomed and reset, the row made current, the rows selected and cleared. A tick alone is not
  evidence: "Select points" and "Deselect points" used to complete on the same selection event.
  Translated from playwright-tests/e2e/tutorials/scatter-plot.test.ts and the TestTrack case
  Apps/tutorials.md (steps 1, 2 and 4: the list, the completion, the hints).

  Serial, because a finished tutorial writes its completion record into the account's settings, which every
  page syncs whole — two features finishing tutorials at once would overwrite each other's records.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "recentViewerSettings" user settings are put back at feature end
    And the "Scatter Plot" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Scatter Plot tutorial
    When user starts the "Scatter Plot" tutorial
    Then the tutorial progress should be 1 of 11
    And "Open scatter plot" tutorial step should be unchecked
    And scatter-plot icon in toolbox should be hinted
    When user clicks on scatter-plot icon in toolbox
    Then the open tableview should have 1 scatter plot viewer
    And "Open scatter plot" tutorial step should be checked
    And the tutorial progress should be 2 of 11

    When user picks "HEIGHT" in the "x" column selector of scatter plot viewer
    Then "Set X to HEIGHT" tutorial step should be checked
    And "xColumnName" property of scatter plot viewer should be "HEIGHT"
    When user picks "WEIGHT" in the "y" column selector of scatter plot viewer
    Then "Set Y to WEIGHT" tutorial step should be checked
    And "yColumnName" property of scatter plot viewer should be "WEIGHT"
    When user picks "AGE" in the "size" column selector of scatter plot viewer
    Then "Set Size to AGE" tutorial step should be checked
    And "sizeColumnName" property of scatter plot viewer should be "AGE"
    When user picks "SEX" in the "color" column selector of scatter plot viewer
    Then "Set Color to SEX" tutorial step should be checked
    And "colorColumnName" property of scatter plot viewer should be "SEX"

    When user remembers the "x axis span" reading of scatter plot viewer
    And user drags a zoom box over the "view" area of scatter plot viewer
    Then "Zoom in" tutorial step should be checked
    And the "x axis span" reading of scatter plot viewer should be lower than remembered
    When user double-clicks on empty plot space of scatter plot viewer
    Then "Double-click to unzoom" tutorial step should be checked
    And the "x axis span" reading of scatter plot viewer should be as remembered

    # 5850 markers sized by AGE overlap: the click makes current the row the plot finds on top under
    # the pointer, so the claim is that row — the plot's own "hovered row" — not the one aimed at
    When user clicks on the "marker of row 11" area of scatter plot viewer
    Then "Click on a point" tutorial step should be checked
    And the "hovered row" reading of scatter plot viewer should be the current row

    When user drags a selection box over the "view" area of scatter plot viewer
    Then "Select points" tutorial step should be checked
    And some rows should be selected
    And "Deselect points" tutorial step should be unchecked
    When user presses Escape in scatter plot viewer
    Then "Deselect points" tutorial step should be checked
    And no rows should be selected

    And the "Scatter Plot" tutorial should be completed
    And the tutorial should have listed 10 steps
    And the tutorial progress should be 11 of 11
    And no hint should be shown
    And no errors should have been logged
