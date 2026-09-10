@journey @viewers @realizes:viewers.line-chart
Feature: Line chart aggregation, whiskers and markers
  What happens when several rows land on the same X position: the chart aggregates them, says which
  function it used, and can draw the spread of what it collapsed as whiskers.
  The old spec proved each of these with an ink delta — "the canvas changed by at least 300 pixels
  within 600 ms", preceded by a 1.5-second poll waiting for the previous change's repaint burst to
  settle. Both are gone: the settle is `isRenderPending` and the claims are `aggregated`,
  `aggregation`, `whiskers`, `markers drawn` and `marker size column`. The gates are the point —
  `whiskers` reads "None" on a chart that is not aggregated no matter what the property says,
  because whiskers need the aggregated frame, and the aggregation-type selector appears only for a
  multi-axis aggregated split.
  Fixture: spgi-100 with X = `Competition assay Date` bucketed by year-quarter — 11 buckets over
  the 80 rows that have a date, so every point is an average of several rows.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | Competition assay Date |
      | xMap         | year quarter           |
      | yColumnNames | Chemical Space X       |
    Then 100 rows should pass the filter
    And the "aggregated" reading of line chart viewer should be "true"
    And the "x categories" reading of line chart viewer should be 11
    And the "aggregation" reading of line chart viewer should be "avg"
    And the "markers drawn" reading of line chart viewer should be 11
    And line chart viewer should report no error

  Scenario: The aggregation function is named, and changing it moves the points
    When user sets "aggrType" property of line chart viewer to "max"
    Then the "aggregation" reading of line chart viewer should be "max"
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be higher than before
    And line chart viewer should have repainted
    When user sets "aggrType" property of line chart viewer to "min"
    Then the "aggregation" reading of line chart viewer should be "min"
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be lower than before
    When user sets "aggrType" property of line chart viewer to "avg"
    Then the "aggregation" reading of line chart viewer should be "avg"
    And no errors should have been logged

  Scenario: Whiskers draw the spread the aggregation collapsed
    Then line chart viewer should not have a "whiskers" area
    And the "whiskers" reading of line chart viewer should be "None"
    When user sets "whiskersType" property of line chart viewer to "Avg | ±StError"
    Then the "whiskers" reading of line chart viewer should be "Avg | ±StError"
    And line chart viewer should have a "whiskers" area
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be higher than before
    And line chart viewer should have repainted
    When user sets "whiskersType" property of line chart viewer to "None"
    Then line chart viewer should not have a "whiskers" area
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be lower than before
    And line chart viewer should have repainted
    And no errors should have been logged

  Scenario: Whiskers on a chart that is not aggregated are not drawn whatever the property says
    When user sets "whiskersType" property of line chart viewer to "Avg | ±StError"
    Then line chart viewer should have a "whiskers" area
    When user sets properties of line chart viewer:
      | xMap        |              |
      | xColumnName | CAST Idea ID |
    Then the "aggregated" reading of line chart viewer should be "false"
    And "whiskersType" property of line chart viewer should be "Avg | ±StError"
    And the "whiskers" reading of line chart viewer should be "None"
    And line chart viewer should not have a "whiskers" area
    When user sets properties of line chart viewer:
      | whiskersType |                        |
      | xColumnName  | Competition assay Date |
      | xMap         | year quarter           |
    Then the "aggregated" reading of line chart viewer should be "true"
    And no errors should have been logged

  Scenario: Markers are drawn unless the chart is told never to draw them
    Then the "markers drawn" reading of line chart viewer should be 11
    When user sets "showMarkers" property of line chart viewer to "Never"
    Then the "markers drawn" reading of line chart viewer should be 0
    And line chart viewer should have less ink than before
    When user sets "showMarkers" property of line chart viewer to "Always"
    Then the "markers drawn" reading of line chart viewer should be 11
    And line chart viewer should have more ink than before
    When user sets "showMarkers" property of line chart viewer to "Auto"
    Then the "markers drawn" reading of line chart viewer should be 11
    And no errors should have been logged

  Scenario: A size column is reported and a marker type change redraws
    Then the "marker size column" reading of line chart viewer should be ""
    When user sets "markerType" property of line chart viewer to "square"
    Then line chart viewer should have repainted
    When user sets "markersSizeColumnName" property of line chart viewer to "Chemical Space Y"
    Then the "marker size column" reading of line chart viewer should be "Chemical Space Y"
    And line chart viewer should have repainted
    And the "markers drawn" reading of line chart viewer should be 11
    When user sets properties of line chart viewer:
      | markersSizeColumnName |        |
      | markerType            | circle |
    Then the "marker size column" reading of line chart viewer should be ""
    And no errors should have been logged

  Scenario: Splitting an aggregated chart aggregates within every category
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then the "lines" reading of line chart viewer should be 5
    And the "categories" reading of line chart viewer should be 5
    And the "markers drawn" reading of line chart viewer should be 32
    And the "x categories" reading of line chart viewer should be 11
    And line chart viewer should have repainted
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "markers drawn" reading of line chart viewer should be 11
    And no errors should have been logged

  Scenario: The aggregation selector shows up only for a multi-axis aggregated split
    Then the "aggr selector shown" reading of line chart viewer should be "false"
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then the "aggr selector shown" reading of line chart viewer should be "false"
    When user sets properties of line chart viewer:
      | yColumnNames | Chemical Space X, TPSA |
      | multiAxis    | true                   |
    Then the "aggr selector shown" reading of line chart viewer should be "true"
    And the "aggregation" reading of line chart viewer should be "avg, avg"
    And the "lines" reading of line chart viewer should be 10
    When user sets "multiAxis" property of line chart viewer to "false"
    Then the "aggr selector shown" reading of line chart viewer should be "false"
    When user sets properties of line chart viewer:
      | yColumnNames     | Chemical Space X |
      | splitColumnNames |                  |
    Then no errors should have been logged
