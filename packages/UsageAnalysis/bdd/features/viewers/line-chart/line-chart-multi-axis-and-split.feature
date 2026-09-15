@journey @viewers @realizes:viewers.line-chart
Feature: Line chart multi-axis layout and splitting into series
  How several Y columns are laid out — one chart box each, or one box with two scales — and how
  many series a split actually draws.
  The whole of the old spec's split half was "the canvas holds more than 38000 painted pixels" with
  a 2-second poll, plus `rowCount === 100` as a "the page is still responsive" probe. The chart now
  says how many series it drew, and the number is the point: `lines` used to be computed as
  `yColumns × categories` with `categories` fixed at 1 unless there was exactly one split column,
  so a two-column split reported 2 where the chart drew 12. On spgi-100 the combinations actually
  present are 5 for Stereo Category, 12 for Stereo × Series (a cartesian product would be 25),
  91 for Stereo × R1, 94 with R2 and 96 with R3 — the numbers this feature asserts.
  Multi-axis is read the same way: `charts` is the boxes laid out, `y axes` the scales drawn, and
  Y Global Scale is visible as the second scale disappearing and the first one stretching to cover
  both columns.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then 100 rows should pass the filter
    And the "charts" reading of line chart viewer should be 1
    And the "lines" reading of line chart viewer should be 1
    And the "categories" reading of line chart viewer should be 1
    And the "split columns" reading of line chart viewer should be 0
    And the "multi axis" reading of line chart viewer should be "false"
    And line chart viewer should report no error

  Scenario: Two Y columns get a chart box each, and Multi Axis folds them into one
    When user sets "yColumnNames" property of line chart viewer to "Chemical Space X, TPSA"
    Then the "charts" reading of line chart viewer should be 2
    And the "lines" reading of line chart viewer should be 2
    And line chart viewer should have a "chart 2" area
    And line chart viewer should have a 'chart "Chemical Space X"' area
    And line chart viewer should have a 'chart "TPSA"' area
    And the "y axes" reading of line chart viewer should be 2
    When user sets "multiAxis" property of line chart viewer to "true"
    Then the "charts" reading of line chart viewer should be 1
    And line chart viewer should not have a "chart 2" area
    And the 'chart "Chemical Space X"' and 'chart "TPSA"' areas of line chart viewer should be the same height
    And the "y axes" reading of line chart viewer should be 2
    And line chart viewer should have a "y2 axis" area
    And line chart viewer should have repainted
    When user sets properties of line chart viewer:
      | multiAxis    | false            |
      | yColumnNames | Chemical Space X |
    Then the "charts" reading of line chart viewer should be 1
    And no errors should have been logged

  Scenario: Y Global Scale replaces the pair of scales with one that covers both columns
    When user sets properties of line chart viewer:
      | yColumnNames | Chemical Space X, TPSA |
      | multiAxis    | true                   |
    Then the "y axes" reading of line chart viewer should be 2
    And line chart viewer should have a "y2 axis" area
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be between 15 and 16
    When user sets "yGlobalScale" property of line chart viewer to "true"
    Then the "y axes" reading of line chart viewer should be 1
    And line chart viewer should not have a "y2 axis" area
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be higher than before
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be between 130 and 133
    And line chart viewer should have repainted
    When user sets "yGlobalScale" property of line chart viewer to "false"
    Then the "y axes" reading of line chart viewer should be 2
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be between 15 and 16
    When user sets properties of line chart viewer:
      | multiAxis    | false            |
      | yColumnNames | Chemical Space X |
    Then no errors should have been logged

  Scenario: One split column draws one series per category
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then the "split columns" reading of line chart viewer should be 1
    And the "lines" reading of line chart viewer should be 5
    And the "categories" reading of line chart viewer should be 5
    And the legend of line chart viewer should list 5 items
    And line chart viewer should have repainted
    And line chart viewer should be painted
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "lines" reading of line chart viewer should be 1
    And the "categories" reading of line chart viewer should be 1
    And no errors should have been logged

  Scenario: A second split column draws the combinations present, not the columns multiplied
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then the "lines" reading of line chart viewer should be 5
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category, Series"
    Then the "split columns" reading of line chart viewer should be 2
    And the "lines" reading of line chart viewer should be 12
    And the "categories" reading of line chart viewer should be 12
    And line chart viewer should be painted
    And line chart viewer should report no error
    When user hovers over the "plot" area of line chart viewer
    Then no errors should have been logged
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "lines" reading of line chart viewer should be 1
    And no errors should have been logged

  Scenario: Every further split column adds only the combinations the rows carry
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category, R1"
    Then the "lines" reading of line chart viewer should be 91
    And the "split columns" reading of line chart viewer should be 2
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category, R1, R2"
    Then the "lines" reading of line chart viewer should be 94
    And the "split columns" reading of line chart viewer should be 3
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category, R1, R2, R3"
    Then the "lines" reading of line chart viewer should be 96
    And the "split columns" reading of line chart viewer should be 4
    And the "rows shown" reading of line chart viewer should be 100
    And line chart viewer should be painted
    And line chart viewer should report no error
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "lines" reading of line chart viewer should be 1
    And the "split columns" reading of line chart viewer should be 0
    And no errors should have been logged

  Scenario: A split multiplies every Y column's series
    When user sets properties of line chart viewer:
      | yColumnNames     | Chemical Space X, TPSA |
      | splitColumnNames | Stereo Category        |
    Then the "charts" reading of line chart viewer should be 2
    And the "lines" reading of line chart viewer should be 10
    And the "categories" reading of line chart viewer should be 5
    When user sets properties of line chart viewer:
      | splitColumnNames |                  |
      | yColumnNames     | Chemical Space X |
    Then the "lines" reading of line chart viewer should be 1
    And no errors should have been logged

  Scenario: Hide other charts on the second chart leaves that column alone
    When user sets "yColumnNames" property of line chart viewer to "Chemical Space X, Chemical Space Y, TPSA"
    Then the "charts" reading of line chart viewer should be 3
    And line chart viewer should have a "chart 3" area
    When user picks "Chemical Space Y > Hide other charts" from the context menu of the "chart 2" area of line chart viewer
    Then the "charts" reading of line chart viewer should be 1
    And the "y columns" reading of line chart viewer should be "Chemical Space Y"
    And line chart viewer should not have a "chart 2" area
    And line chart viewer should be painted
    When user sets "yColumnNames" property of line chart viewer to "Chemical Space X"
    Then the "y columns" reading of line chart viewer should be "Chemical Space X"
    And no errors should have been logged
