@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: a bar chart set up from the chart, its context menu and its properties. The
  category is picked in the selector beside the bars (DIS_POP: six bars); the value and its aggregation
  in the one under them (avg of WEIGHT); stacking needs an additive aggregation, so the aggregation goes
  back to count before SEX is picked in the Stack selector at the top (twelve segments). The bars come
  sorted by value, largest first (RA above AS); the context menu sorts them by name instead (Order > by
  category, then asc: AS above RA) and turns them upright (Orientation > vertical: AS left of RA). The
  gear's Y Axis > Relative Values shows each bar as 100%. The aggregation list is a native select,
  picked by the step in howto-viewers.ts. The chart is added on RACE and given most of the view's width
  off camera; demog is opened off camera; nothing is made on the server.

  Scenario: Configure a bar chart
    Pick the category, the value and how to aggregate it, stack the bars, then sort them and turn them upright.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a bar chart viewer with:
      | Split | RACE |
    And bar chart viewer takes 65% of the space it shares with the grid
    When the video says "A bar chart compares one number across categories"
    # caption: The categories: the selector beside the bars
    And user picks "DIS_POP" in the "split" column selector of bar chart viewer
    Then the "bars" reading of bar chart viewer should be 6
    # caption: The number: how to aggregate, and which column
    When user picks the "avg" aggregation in the "value" column selector of bar chart viewer
    And user picks "WEIGHT" in the "value" column selector of bar chart viewer
    Then "Value Aggr Type" property of bar chart viewer should be "avg"
    And "Value" property of bar chart viewer should be "WEIGHT"
    # caption: Stacking needs a count or a sum
    When user picks the "count" aggregation in the "value" column selector of bar chart viewer
    # caption: Stack: split each bar by another category
    And user picks "SEX" in the "stack" column selector of bar chart viewer
    Then the "stack segments" reading of bar chart viewer should be 12
    And the "bar RA" area of bar chart viewer should lie above the "bar AS" area
    # caption: Right-click to sort the bars by name
    When user picks "Order > by category" from the context menu of the "view" area of bar chart viewer
    And user picks "Order > asc" from the open menu
    And user closes the context menu
    Then the "bar AS" area of bar chart viewer should lie above the "bar RA" area
    # caption: or to turn them upright
    When user picks "Orientation > vertical" from the context menu of the "view" area of bar chart viewer
    And user closes the context menu
    Then the "bar AS" area of bar chart viewer should lie to the left of the "bar RA" area
    # caption: The gear opens every setting, grouped in sections
    When user hovers over the "x axis" area of bar chart viewer
    And user clicks on settings icon of bar chart viewer
    Then the context panel should show "Bar chart"
    # caption: Y Axis: show each bar as 100%
    When user takes a snapshot of bar chart viewer
    And user expands "Y Axis" category in context panel
    And user checks "Relative Values" property in context panel
    Then "Relative Values" property of bar chart viewer should be "true"
    And bar chart viewer should have repainted
    And no errors should have been logged
