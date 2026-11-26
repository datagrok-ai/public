@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: a line chart set up from the chart, its properties and its context menu. It starts
  on a number (AGE) as X; STARTED, a date, is picked in the X selector, and the date list beside it
  groups the dates by month (the chart then aggregates each month's rows); the split selector at the top
  left adds one line per SEX. The gear then adds HEIGHT as a second Y column in the Y Axis section's
  column list (a chart each), and Data > Multi Axis puts both into one chart, which then shows the "Y
  Aggr Type" list above the plot: the median there aggregates every line. (The context menu's Data >
  Aggregation is not used: it changes the lines but leaves that list blank and Aggr Type at avg, so a
  column added later is averaged.) Last, Lines > Interpolation draws smooth lines. The date grouping and
  the Y Aggr Type list are native selects, picked by the steps in howto-viewers.ts. The chart is added
  (X AGE, Y WEIGHT) under the grid and given most of the view's height off camera; demog is opened off
  camera; nothing is made on the server.

  Scenario: Configure a line chart
    Put a date on X and group it by month, split by a category, add a second Y column, change the aggregation.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | WEIGHT |
    And line chart viewer takes 65% of the space it shares with the grid
    When the video says "A line chart shows how a value changes along X"
    # caption: X can be a number or a date
    And user picks "STARTED" in the "x" column selector of line chart viewer
    Then the "x column" reading of line chart viewer should be "STARTED"
    # caption: Group the dates, here by month
    When user picks the "year month" aggregation in the "x" column selector of line chart viewer
    Then the "aggregated" reading of line chart viewer should be "true"
    And the "aggregation" reading of line chart viewer should be "avg"
    # caption: Split: one line for each category
    When user adds "SEX" to the splits of line chart viewer
    Then the "lines" reading of line chart viewer should be 2
    # caption: The gear opens every setting, grouped in sections
    When user hovers over the "x axis" area of line chart viewer
    And user clicks on settings icon of line chart viewer
    Then the context panel should show "Line chart"
    # caption: Y Axis: add more columns, each gets its own chart
    When user expands "Y Axis" category in context panel
    And user clicks on "..." button in "Y" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "HEIGHT" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "charts" reading of line chart viewer should be 2
    # caption: Data: Multi Axis draws them in one chart
    When user checks "Multi Axis" property in context panel
    Then the "charts" reading of line chart viewer should be 1
    And the "y axes" reading of line chart viewer should be 2
    # caption: Y Aggr Type: how each month's rows are combined
    When user picks "med" in the Y Aggr Type selector of line chart viewer
    Then the "aggregation" reading of line chart viewer should be "med, med"
    And "Aggr Type" property of line chart viewer should be "med"
    # caption: Lines: smooth the lines
    When user expands "Lines" category in context panel
    And user selects "Spline" in "Interpolation" property in context panel
    Then "Interpolation" property of line chart viewer should be "Spline"
    And no errors should have been logged
