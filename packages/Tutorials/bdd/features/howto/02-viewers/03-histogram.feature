@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: a histogram set up from the plot, then from its properties. The value column is
  picked in the selector under the plot (AGE); a hovered bar shows its count; the bins slider at the
  top of the plot (shown under the pointer) changes the number of bins; the range handles under the
  bars filter the table, and a double-click on the range slider releases it; the split selector at the
  top right draws one line per category (SEX). Then the gear opens the properties: Category > Split
  Stack brings the bars back as one segment per category, X Axis > Show Values writes the counts over
  the bars. The histogram is added on WEIGHT and given most of the view's width off camera; demog is
  opened off camera; nothing is made on the server.

  Scenario: Configure a histogram
    Pick the column, change the bins, filter with the range slider, split by a category, then use the settings.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a histogram viewer with:
      | Value | WEIGHT |
    And histogram viewer takes 65% of the space it shares with the grid
    When the video says "A histogram shows how the values of one column are spread"
    # caption: Pick the column under the plot
    And user picks "AGE" in the "value" column selector of histogram viewer
    Then "Value" property of histogram viewer should be "AGE"
    # caption: The slider at the top sets the number of bins
    When user drags the slider of histogram viewer to 12
    Then the "bins shown" reading of histogram viewer should be between 10 and 14
    # caption: Hover a bar to see how many rows it holds
    When user hovers over the "bin 6" area of histogram viewer
    # caption: The handles under the bars filter the table
    When user hovers over the "view" area of histogram viewer
    And user drags the "range max handle" area of histogram viewer by 400 pixels to the left
    Then fewer than 5850 rows should pass the filter
    # caption: Double-click the slider to show all rows again
    When user double-clicks on the "range slider" area of histogram viewer
    Then all rows should pass the filter
    # caption: Split by a category: one line for each
    When user picks "SEX" in the "split" column selector of histogram viewer
    Then histogram viewer should have a "line F" area
    And histogram viewer should have a "line M" area
    # caption: The gear opens every setting, grouped in sections
    When user hovers over the "x axis" area of histogram viewer
    And user clicks on settings icon of histogram viewer
    Then the context panel should show "Histogram"
    # caption: Category: stack the categories in each bar
    When user expands "Category" category in context panel
    And user checks "Split Stack" property in context panel
    Then histogram viewer should have a "bin 8 | F" area
    # caption: X Axis: write the count over each bar
    When user expands "X Axis" category in context panel
    And user checks "Show Values" property in context panel
    Then histogram viewer should have a "bin labels" area
    And no errors should have been logged
