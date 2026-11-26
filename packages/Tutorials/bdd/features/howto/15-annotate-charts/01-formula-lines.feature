@guide @tutorial
Feature: How-to · Annotations
  For the maintainer: formula lines on a scatter plot of demog-1000, HEIGHT (cm) against WEIGHT (kg); 872 rows have
  both. The Y axis's context menu, Annotations > Add Line, adds a horizontal line at the median WEIGHT and opens the
  PowerPack Formula Lines dialog on it; its Value is set to 100 (the line takes the title "WEIGHT = 100"). Add new >
  Line then takes a formula of both columns: the weight at a BMI of 25, ${WEIGHT} = 25 * (${HEIGHT} / 100) ^ 2,
  titled "BMI 25" and colored red in the Format pane. Every point above that curve has a BMI over 25 (544 of the 872). OK keeps both on the viewer (its
  Formula Lines property, in the browser only); the plot's Tools > Formula Lines... opens the same dialog listing
  them, and CANCEL closes it. demog-1000 is opened and the plot added off camera; nothing is made on the server.

  Scenario: Draw reference lines with formula lines
    A formula line marks a value or a curve on a chart: a limit, a target, a known relationship.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X | HEIGHT |
      | Y | WEIGHT |
    And scatter plot viewer takes 65% of the space it shares with the grid
    When the video says "Formula lines mark a value or a curve on a chart"
    # caption: Right-click an axis: Annotations > Add Line
    And user picks "Annotations > Add Line" from the context menu of the "y axis" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    # caption: Type the value the line marks
    When user enters "100" into "Value" text input in "Formula Lines" dialog
    # caption: Add new > Line takes a formula of X and Y
    And user clicks on "Add new" button in "Formula Lines" dialog
    And user picks "Line" from the open menu
    # caption: Here: the weight at a BMI of 25
    And user enters "${WEIGHT} = 25 * (${HEIGHT} / 100) ^ 2" into text area in "Line" section
    # caption: The title is written on the line
    And user enters "BMI 25" into "Title" text input in "Formula Lines" dialog
    # caption: Format: the line's color, style and width
    And user enters "#d62728" into "Color" input in "Formula Lines" dialog
    And user clicks on OK button in "Formula Lines" dialog
    # caption: Points above the curve have a BMI over 25
    Then scatter plot viewer should draw 2 formula lines
    And the "view" area of scatter plot viewer should contain the color "#d62728"
    # caption: Tools > Formula Lines... lists and edits them all
    When user picks "Tools > Formula Lines..." from the context menu of the "view" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on CANCEL button in "Formula Lines" dialog
    Then the formula lines of scatter plot viewer should hold "${WEIGHT} = 100" over the ranges ".."
    And the formula lines of scatter plot viewer should hold "${WEIGHT} = 25 * (${HEIGHT} / 100) ^ 2" over the ranges ".."
    And no errors should have been logged
