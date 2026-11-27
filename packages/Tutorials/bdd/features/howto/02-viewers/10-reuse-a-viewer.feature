@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: two ways to reuse a viewer you have set up. Two scatter plots of demog-1000: the
  first (HEIGHT against WEIGHT) has hollow markers and a histogram along each axis, the second (AGE
  against HEIGHT) is as it comes. Pick Up / Apply > Pick Up on the first remembers its settings; Pick
  Up / Apply > Apply Style Settings on the second copies the look (hollow markers, both histograms) and
  keeps its own columns (Apply would copy the columns too). General > Clone then adds a copy of the
  second plot next to it. Nothing is saved as a default (Set as Default is not used); the plots are
  added and styled and demog-1000 opened off camera; nothing is made on the server.

  Scenario: Copy a viewer or its look
    Copy one viewer's look onto another with Pick Up and Apply, or clone the viewer.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X                | HEIGHT        |
      | Y                | WEIGHT        |
      | Marker Type      | circle border |
      | Show X Histogram | true          |
      | Show Y Histogram | true          |
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | HEIGHT |
    Then the open tableview should have 2 scatter plot viewers
    And second scatter plot viewer should not have an "x histogram" area
    When the video says "Reuse the look of a viewer you have set up"
    # caption: Right-click > Pick Up / Apply > Pick Up remembers its settings
    And user picks "Pick Up / Apply > Pick Up" from the context menu of the "view" area of first scatter plot viewer
    # caption: Apply Style Settings copies the look, not the columns
    And user picks "Pick Up / Apply > Apply Style Settings" from the context menu of the "view" area of second scatter plot viewer
    Then second scatter plot viewer should have an "x histogram" area
    And second scatter plot viewer should have a "y histogram" area
    And "Marker Type" property of second scatter plot viewer should be "circle border"
    And "X" property of second scatter plot viewer should be "AGE"
    # caption: General > Clone adds a copy of a viewer
    When user picks "General > Clone" from the context menu of the "view" area of second scatter plot viewer
    Then the open tableview should have 3 scatter plot viewers
    And "X" property of third scatter plot viewer should be "AGE"
    And third scatter plot viewer should have an "x histogram" area
    And no errors should have been logged
