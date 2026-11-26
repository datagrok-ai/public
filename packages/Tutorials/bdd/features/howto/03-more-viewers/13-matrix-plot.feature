@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The matrix plot from the toolbox on iris: a cell per pair
  of numeric columns with a histogram on the diagonal, then the settings icon: scatter plots in the cells instead
  of density plots (Misc), the row-number column "col 1" taken off both axes in the column dialogs (Data), and
  the cells' scatter plots colored by Species from the inner viewer's tab. iris is opened off camera; nothing is
  made on the server.

  Scenario: See every pair of columns with a matrix plot
    A matrix plot draws a small plot for every pair of numeric columns, so you can scan all relationships at once.

    Given user is logged in
    And simple mode is off
    And user opens iris dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A matrix plot: every pair of numeric columns in one grid"
    And user clicks on matrix plot icon on toolbox
    Then matrix plot viewer should be added to the open tableview
    # caption: Hide the toolbox to give the plot room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    When the video says "Each cell plots one column against another; the diagonal shows histograms" at the top
    And user hovers over the "cell body Petal.Length x Sepal.Length" area of matrix plot viewer
    # caption: The rest is in the settings
    And user hovers over matrix plot viewer
    And user clicks on settings icon of matrix plot viewer
    Then the context panel should show "Matrix plot"
    # caption: Misc: scatter plots instead of density plots
    When user expands "Misc" category in context panel
    And user selects "Scatter plot" in "Cell Plot Type" property in context panel
    Then the "cell viewer type" reading of matrix plot viewer should be "Scatter plot"
    # caption: Data: choose the columns for each axis
    When user clicks on "..." button in "X" property in context panel
    And user toggles the "col 1" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user clicks on "..." button in "Y" property in context panel
    And user toggles the "col 1" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "cells" reading of matrix plot viewer should be 16
    # caption: The second tab sets up the plot inside every cell
    When user clicks on "Scatter plot" tab in context panel
    And user expands "Color" category in context panel
    And user selects "Species" in "Color" property in context panel
    Then "colorColumnName" inner property of matrix plot viewer should be "Species"
    And no errors should have been logged
