@guide @help:visualize/viewers
Feature: Build a dashboard on molecular data
  A guide: the answer to "how do I build a dashboard on molecular data, with several viewers set up
  the way I want?". Every viewer comes from the Toolbox's Viewers pane, docks next to the grid, and
  is set up where a person sets it up: the column selectors drawn on the viewer itself (X, Y and
  Color on a scatter plot, the value of a histogram, the split of a bar chart, the category of a
  pie chart), the settings icon in its corner for everything else, and the ribbon's filter icon for
  the filter panel. Demo: spgi-100, a hundred structures with their chemical-space coordinates,
  series and properties.

  Scenario: Add viewers from the toolbox and set each one up
    Given user is logged in
    And no project named "Molecular dashboard" is on the server
    And user opens spgi dataset
    When user opens toolbox
    And user clicks on scatter plot icon on toolbox
    Then scatter plot viewer should be added to the open tableview
    When user picks "Chemical Space X" in the "x" column selector of scatter plot viewer
    And user picks "Chemical Space Y" in the "y" column selector of scatter plot viewer
    And user picks "Series" in the "color" column selector of scatter plot viewer
    Then scatter plot viewer should be painted in at least 3 colors
    When user clicks on histogram icon on toolbox
    And user selects "Average Mass" in Value column input in histogram viewer
    And user clicks on bar chart icon on toolbox
    And user picks "Series" in the "split" column selector of bar chart viewer
    And user clicks on pie chart icon on toolbox
    And user picks "Stereo Category" in the "category" column selector of pie chart viewer
    And user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user clicks on settings icon of bar chart viewer
    Given "Y Axis" category in context panel is expanded
    When user selects "avg" in "Value Aggr Type" property in context panel
    And user selects "Average Mass" in "Value" property in context panel
    Then the current view should hold at least 5 viewers
    When user clicks on Save button in toolbar
    And user enters "Molecular dashboard" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then 1 project named "Molecular dashboard" should be on the server
    When user clicks on CANCEL button in "Share Molecular dashboard" dialog
    Then no errors should have been logged
