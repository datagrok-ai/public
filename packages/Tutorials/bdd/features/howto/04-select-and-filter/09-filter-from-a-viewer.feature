@guide @tutorial
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Viewers that filter: On Click > Filter in a bar chart's
  context menu makes a click on a bar filter the table to that category (and Row Source All keeps every bar on
  the chart, so another bar can be clicked), a click on empty space lets go, and on a scatter plot, whose Zoom
  and Filter is "filter by zoom" by default, dragging the end of an axis slider zooms in and filters the table to
  the points left in view, until Reset View. demog (5850 rows: DIS_POP RA 2550, UC 728) with a bar chart split by
  DIS_POP and a scatter plot of AGE and WEIGHT is opened off camera; nothing is made on the server.

  Scenario: Filter from a viewer
    Click a bar to filter by its category, or zoom a scatter plot to filter by what stays in view.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And user adds a bar chart viewer with:
      | Split | DIS_POP |
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    Then all rows should pass the filter
    # caption: A bar chart can filter on click instead of selecting
    When user picks "On Click > Filter" from the context menu of bar chart viewer
    And user closes the context menu
    Then "On Click" property of bar chart viewer should be "Filter"
    # caption: Click a bar to keep only its rows
    When user clicks on the "bar RA" area of bar chart viewer
    Then 2550 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" is "RA"
    When user clicks on the "bar UC" area of bar chart viewer
    Then 728 rows should pass the filter
    # caption: Click empty space to let go
    When user clicks on empty plot space of bar chart viewer
    Then all rows should pass the filter
    # caption: Drag the end of an axis slider to zoom in
    When user drags the max handle of the "x" range slider of scatter plot viewer by -350 pixels
    Then fewer than 4000 rows should pass the filter
    When user moves the pointer away from scatter plot viewer
    And the video says "Zooming a scatter plot filters the table to the points in view"
    # caption: Reset View brings every row back
    And user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And no errors should have been logged
