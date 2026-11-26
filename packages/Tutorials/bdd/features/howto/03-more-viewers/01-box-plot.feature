@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The box plot from the toolbox (the toolbox then hidden on
  camera with its sidebar icon, so the plot gets the room without a jump), its Value and Category
  picked on the plot, the p-value under it, then the settings icon: a second category (Data), the statistics
  under the plot (Statistics), the violin style (Style) and the box color by a numeric column (Color).
  demog is opened off camera; nothing is made on the server.

  Scenario: Compare groups with a box plot
    A box plot shows how a number spreads in each group: the median, the middle half and the outliers.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A box plot compares a number across groups"
    And user clicks on box plot icon on toolbox
    Then box plot viewer should be added to the open tableview
    # caption: Hide the toolbox to give the plot room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Value: the number to compare
    When user picks "HEIGHT" in the "Value" column selector of box plot viewer
    # caption: Category: the groups
    And user picks "RACE" in the "Category 1" column selector of box plot viewer
    When the video says "Under the plot: min, max, average and median of each group" at the top
    And user hovers over the "stats" area of box plot viewer
    # caption: The p-value tells whether the groups really differ
    And user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "Alexander and Govern"
    # caption: Everything else is in the settings
    When user hovers over box plot viewer
    And user clicks on settings icon of box plot viewer
    Then the context panel should show "Box plot"
    # caption: Data: split each group again by a second category
    When user selects "SEX" in "Category2" property in context panel
    Then box plot viewer should have a "category Asian | F" area
    # caption: Statistics: choose what the table under the plot shows
    When user expands "Statistics" category in context panel
    And user checks "Show Q1" property in context panel
    And user checks "Show Q3" property in context panel
    Then box plot viewer should have a "stats" area
    When user collapses "Statistics" category in context panel
    # caption: Color: color the dots by any column
    When user expands "Color" category in context panel
    And user selects "AGE" in "Marker Color" property in context panel
    Then box plot viewer should have a "color scale" area
    # caption: Style: a violin shows the shape of each distribution
    When user expands "Style" category in context panel
    And user selects "violin" in "Plot Style" property in context panel
    Then "Plot Style" property of box plot viewer should be "violin"
    And no errors should have been logged
