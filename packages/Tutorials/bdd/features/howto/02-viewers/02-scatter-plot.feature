@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: a scatter plot set up from the plot itself, then from its properties, on cars (30
  rows, so the labels can be read). X and Y are picked in the selectors on the axes (horsepower against
  price), Color and Size in the ones at the top right, which show under the pointer (city.mpg and
  curb.weight). The wheel zooms, and since the plot's Zoom and Filter is "filter by zoom" by default,
  the table keeps only the rows in view; a plain drag pans and a double-click on empty space brings
  everything back. The context menu's Tools > Show Regression Line draws the fit. The gear then opens
  the properties in the context panel, and one setting is changed in each of four sections: Marker (the
  marker type), Labels (the model names, picked in the column dialog), X Axis (a histogram along the
  axis) and Statistics (Pearson's r beside the fit). The plot is added off camera on wheel.base and
  length, so neither pick plots a column against itself; cars is opened off camera; nothing is made on
  the server.

  Scenario: Configure a scatter plot
    Pick the columns on the plot, zoom and pan, add a regression line, then change its settings from the gear icon.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X | wheel.base |
      | Y | length     |
    And scatter plot viewer takes 65% of the space it shares with the grid
    When the video says "A scatter plot draws one dot per row"
    # caption: Pick X and Y right on the axes
    And user picks "horsepower" in the "x" column selector of scatter plot viewer
    And user picks "price" in the "y" column selector of scatter plot viewer
    Then "X" property of scatter plot viewer should be "horsepower"
    And "Y" property of scatter plot viewer should be "price"
    # caption: Color and Size: the selectors at the top right
    When user picks "city.mpg" in the "color" column selector of scatter plot viewer
    Then scatter plot viewer should have a "color scale" area
    When user picks "curb.weight" in the "size" column selector of scatter plot viewer
    Then "Size" property of scatter plot viewer should be "curb.weight"
    # caption: The wheel zooms, and the table keeps the rows in view
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 30 rows should pass the filter
    # caption: Drag to move around
    When user drags across the "view" area of scatter plot viewer
    # caption: Double-click an empty spot to see all the points again
    And user double-clicks on empty plot space of scatter plot viewer
    Then 30 rows should pass the filter
    # caption: Right-click for more, like a regression line
    When user picks "Tools > Show Regression Line" from the context menu of the "view" area of scatter plot viewer
    And user closes the context menu
    Then scatter plot viewer should have a "regression line" area
    # caption: The gear opens every setting, grouped in sections
    When user hovers over empty plot space of scatter plot viewer
    And user clicks on settings icon of scatter plot viewer
    Then the context panel should show "Scatter plot"
    # caption: Marker: the shape of the dots
    When user expands "Marker" category in context panel
    And user selects "square" in "Marker Type" property in context panel
    Then "Marker Type" property of scatter plot viewer should be "square"
    # caption: Labels: write a column's value next to each dot
    When user expands "Labels" category in context panel
    And user clicks on "..." button in "Label Columns" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "model" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Label Columns" property of scatter plot viewer should contain "model"
    And the "labels shown" reading of scatter plot viewer should be at least 10
    # caption: X Axis: a histogram along the axis
    When user expands "X Axis" category in context panel
    And user checks "Show X Histogram" property in context panel
    Then scatter plot viewer should have an "x histogram" area
    # caption: Statistics: what to show beside the fit
    When user takes a snapshot of scatter plot viewer
    And user expands "Statistics" category in context panel
    And user checks "Show Pearson Correlation" property in context panel
    Then the "regression stats" area of scatter plot viewer should be wider than before
    And no errors should have been logged
