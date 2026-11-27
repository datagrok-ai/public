@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The correlation plot from the toolbox on cars: a cell
  hovered (its scatter plot in the tooltip), then the settings icon: Spearman instead of Pearson (Misc), the Y
  axis cut down to city.mpg and price in the column dialog (Data; two columns, since a single Y column drew its
  cells without colors on 2026-11-22), and a cell double-clicked into a real scatter plot of the pair.
  cars is opened off camera; nothing is made on the server.

  Scenario: Find related columns with a correlation plot
    A correlation plot shows how strongly every pair of numeric columns moves together.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A correlation plot compares every pair of numeric columns"
    And user clicks on correlation plot icon on toolbox
    Then correlation plot viewer should be added to the open tableview
    # caption: Hide the toolbox to give the plot room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    When the video says "Red: the two rise together. Blue: one falls as the other rises" at the bottom
    # caption: Hover a cell to see that pair as a scatter plot
    And user hovers over the "cell horsepower x price" area of correlation plot viewer
    # caption: The rest is in the settings
    When user hovers over correlation plot viewer
    And user clicks on settings icon of correlation plot viewer
    Then the context panel should show "Correlation plot"
    # caption: Misc: Spearman compares ranks, so outliers count less
    When user expands "Misc" category in context panel
    And user selects "Spearman" in "Correlation Type" property in context panel
    Then the "correlation type" reading of correlation plot viewer should be "Spearman"
    # caption: Data: keep just mileage and price on one axis
    When user clicks on "..." button in "Y" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    And user toggles the "city.mpg" column in the column list of "Select columns..." dialog
    And user toggles the "price" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "y columns" reading of correlation plot viewer should be "city.mpg, price"
    # caption: Double-click a cell to open the pair as a scatter plot
    When user double-clicks on the "cell horsepower x price" area of correlation plot viewer
    Then the open tableview should have 1 scatter plot viewer
    When user moves the pointer away from correlation plot viewer
    Then no errors should have been logged
