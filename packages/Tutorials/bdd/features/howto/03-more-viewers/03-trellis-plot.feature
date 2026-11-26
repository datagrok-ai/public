@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The trellis plot from the toolbox, split by SEX across
  and RACE down with the column selectors on the plot, a scatter plot picked as the chart in every cell, then the
  settings icon. demog is opened off camera; nothing is made on the server.

  Scenario: Small multiples with a trellis plot
    A trellis plot draws the same chart once for every category, so the groups sit side by side.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A trellis plot draws one small chart per category"
    And user clicks on trellis plot icon on toolbox
    Then trellis plot viewer should be added to the open tableview
    # caption: Hide the toolbox to give the plot room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: One column splits the cells across
    When user picks "SEX" in the column selector at the "x selector 1" area of trellis plot viewer
    # caption: Another splits them down
    And user picks "RACE" in the column selector at the "y selector 1" area of trellis plot viewer
    Then the "cells" reading of trellis plot viewer should be 8
    # caption: Choose the chart that goes in each cell
    When user picks "Scatter plot" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"
    # caption: Its columns apply to every cell
    When user picks "AGE" in the "color" column selector of trellis plot viewer
    # caption: The rest is in the settings
    When user hovers over trellis plot viewer
    And user clicks on settings icon of trellis plot viewer
    Then the context panel should show "Trellis plot"
    # caption: Axes: one scale for every cell, so they compare
    When user expands "Axes" category in context panel
    And user checks "Global Scale" property in context panel
    # caption: Misc: click a cell to filter the table to it
    When user expands "Misc" category in context panel
    And user selects "Filter" in "On Click" property in context panel
    And user clicks on the "cell M | Black" area of trellis plot viewer
    Then 53 rows should pass the filter
    # caption: The second tab holds the settings of the chart inside
    When user clicks on "Scatter plot" tab in context panel
    And user expands "Statistics" category in context panel
    And user checks "Show Regression Line" property in context panel
    Then "showRegressionLine" inner property of trellis plot viewer should be "true"
    And no errors should have been logged
