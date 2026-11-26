@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: the two ways to add a viewer to a table, then arranging them. The toolbox's
  Viewers icons add a scatter plot and a histogram with one click each; the toolbar's Add viewer icon
  opens the gallery of every viewer, where a tag (Trend) narrows the cards to the line chart, and a
  search ("bar") finds the bar chart. Then the viewers are arranged by their title bars: the bar chart
  dragged onto the bottom edge of the view runs under the grid and the other viewers, the scatter plot
  goes full screen from its title bar's arrows icon and comes back with Alt+F (Esc in a real browser;
  the filming browser ignores it), and the histogram is closed with the x on its title bar. The title
  bar's icons show only under the pointer, hence the hovers. demog is opened off camera; nothing is
  made on the server.

  Scenario: Add viewers to a table
    Add a viewer from the toolbox or from the gallery of all viewers, then arrange the viewers around the grid.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Viewers are charts of the open table"
    # caption: The toolbox: the common viewers, one click each
    And user clicks on toolbox tab
    Then toolbox should be visible
    When user clicks on "scatter plot" icon in toolbox
    Then scatter plot viewer should be visible
    When user clicks on "histogram" icon in toolbox
    Then histogram viewer should be visible
    When user clicks on toolbox tab
    # caption: Add viewer opens the gallery of every viewer
    And user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    # caption: Tags group the viewers by what they show
    When user clicks on "Trend" viewer tag
    And user clicks on "Line chart" viewer card
    Then "Add Viewer" dialog should be absent
    And line chart viewer should be visible
    When user clicks on "Add viewer" icon in toolbar
    # caption: Or search by name
    And user types "bar" into viewer gallery search
    And user clicks on "Bar chart" viewer card
    Then bar chart viewer should be visible
    # caption: Drag a viewer by its title bar to dock it elsewhere
    When user docks bar chart viewer to the bottom edge of the view
    Then bar chart viewer should be docked along the bottom edge of the view
    # caption: The arrows on the title bar: full screen
    When user hovers over empty plot space of scatter plot viewer
    And user clicks on "expand-arrows" icon in "Scatter plot" dock panel
    # caption: Alt+F or Esc goes back
    And user presses Alt+F in scatter plot viewer
    # caption: The x on the title bar closes a viewer
    And user hovers over empty plot space of histogram viewer
    And user clicks on close icon of histogram viewer
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 1 scatter plot viewer
    And bar chart viewer should be docked along the bottom edge of the view
    And no errors should have been logged
