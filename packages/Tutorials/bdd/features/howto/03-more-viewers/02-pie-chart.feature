@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The pie chart from the toolbox, its category picked on the
  chart, a slice hovered and clicked (its 650 rows selected), then the settings icon: the largest slice first
  (Data), the donut (Style), the row count on each slice (Misc).
  demog is opened off camera; nothing is made on the server.

  Scenario: Show shares with a pie chart
    A pie chart splits the rows by one category; each slice is that category's share.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A pie chart shows each category's share of the rows"
    And user clicks on pie chart icon on toolbox
    Then pie chart viewer should be added to the open tableview
    # caption: Hide the toolbox to give the chart room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Pick the category right on the chart
    When user picks "SEVERITY" in the "category" column selector of pie chart viewer
    Then the "slices" reading of pie chart viewer should be 5
    # caption: Hover a slice for its row count
    When user hovers over the "slice High" area of pie chart viewer
    # caption: Click a slice to select its rows
    And user clicks on the "slice High" area of pie chart viewer
    Then 650 rows should be selected
    # caption: The rest is in the settings
    When user hovers over pie chart viewer
    And user clicks on settings icon of pie chart viewer
    Then the context panel should show "Pie chart"
    # caption: Data: sort the slices, largest first
    When user selects "desc" in "Pie Sort Order" property in context panel
    Then the "start angle of None" reading of pie chart viewer should be 0
    # caption: Style: turn the pie into a donut
    When user expands "Style" category in context panel
    And user selects "Donut" in "Mode" property in context panel
    Then "Mode" property of pie chart viewer should be "Donut"
    # caption: Misc: show the row count on each slice
    When user expands "Misc" category in context panel
    And user checks "Show Value" property in context panel
    Then no errors should have been logged
