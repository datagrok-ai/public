@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The PC plot (parallel coordinates) from the toolbox on
  iris: colored by Species on the plot, a line hovered, the Petal.Length axis narrowed by dragging its top handle
  (the table filtered to the rows in range), then the settings icon: the filtered-out lines drawn (Data), the
  row-number column "col 1" taken off the axes and one scale for the four measurements in cm (Y Axis), and
  Reset View from the plot's menu. iris is opened off camera; nothing is made on the server.

  Scenario: Compare many columns with a PC plot
    A PC plot draws every row as a line across one axis per column, so patterns across many columns show up.

    Given user is logged in
    And simple mode is off
    And user opens iris dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A PC plot draws each row as a line across many columns"
    And user clicks on pc plot icon on toolbox
    Then pc plot viewer should be added to the open tableview
    # caption: Hide the toolbox to give the plot room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Color the lines by a column
    When user picks "Species" in the "color" column selector of pc plot viewer
    # caption: Drag an axis handle to keep only the rows in range
    And user hovers over pc plot viewer
    And user drags the max handle of the "Petal.Length" range slider of pc plot viewer by 300 pixels
    Then fewer than 150 rows should pass the filter
    # caption: The rest is in the settings
    When user clicks on settings icon of pc plot viewer
    Then the context panel should show "PC Plot"
    # caption: Data: keep the filtered-out rows in view, in grey
    When user checks "Show Filtered Out Lines" property in context panel
    Then the "filtered out lines drawn" reading of pc plot viewer should be at least 1
    # caption: Y Axis: choose the columns to draw
    When user expands "Y Axis" category in context panel
    And user clicks on "..." button in "Column Names" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "col 1" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "axis order" reading of pc plot viewer should be "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width"
    # caption: All four are in cm, so one scale fits them all
    When user unchecks "Normalize Each Column" property in context panel
    Then the "normalization" reading of pc plot viewer should be "global"
    # caption: Reset View clears the axis filters
    When user picks "Reset View" from the context menu of pc plot viewer
    Then all rows should pass the filter
    And no errors should have been logged
