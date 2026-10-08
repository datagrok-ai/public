@journey @viewers @realizes:viewers.tooltips
Feature: The line chart's aggregated tooltip with a split
  A line chart whose X has repeated values aggregates its points, and its Tooltip > Edit... opens
  Edit Aggregated Tooltip, where each row adds an aggregation of a column to the tooltip. With two
  aggregations configured and the chart split by a column into several lines, hovering a point
  shows the configured aggregations (GitHub #2571). Translated from the TestTrack case
  Tooltips/line-chart-aggregated-tooltip, on spgi-100 (the first 100 rows of SPGI) instead of the
  whole SPGI: X = Chemist 521, Y = CAST Idea ID, the aggregations concat unique(Stereo Category) and
  min(Average Mass), split by Stereo Category (five lines: four categories and the empty one). The
  pointer rests on the first point the split chart reports (`point of row <r>`, a row of the
  aggregated frame), whose tooltip lists the two aggregations.

  The case's title also names an in-viewer filter, which none of its steps sets; not translated.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And user opens spgi dataset
    And user adds a line chart viewer
    And user sets properties of line chart viewer:
      | X | Chemist 521  |
      | Y | CAST Idea ID |

  Scenario: Edit Aggregated Tooltip adds two aggregations
    Then the "aggregated" reading of line chart viewer should be "true"
    When user picks "Tooltip > Edit..." from the context menu of line chart viewer
    Then "Edit Aggregated Tooltip" dialog should be visible
    When user clicks on first button in "Edit Aggregated Tooltip" dialog
    And user picks "Stereo Category" in the column selector of the "Edit Aggregated Tooltip" dialog
    And user selects "concat unique" in choice input in "Edit Aggregated Tooltip" dialog
    And user clicks on "Add aggregation" button in "Edit Aggregated Tooltip" dialog
    And user selects "Average Mass" in second column selector in "Edit Aggregated Tooltip" dialog
    And user selects "min" in second choice input in "Edit Aggregated Tooltip" dialog
    And user clicks on OK button in "Edit Aggregated Tooltip" dialog
    Then the "Edit Aggregated Tooltip" dialog should close
    And "aggTooltipColumns" property of line chart viewer should be "concat unique(Stereo Category)\nmin(Average Mass)"
    And no errors should have been logged

  Scenario: Split by Stereo Category, a point's tooltip shows the two aggregations
    When user sets "Split" property of line chart viewer to "Stereo Category"
    Then the "split columns" reading of line chart viewer should be 1
    And the "lines" reading of line chart viewer should be 5
    And the "markers drawn" reading of line chart viewer should be 32
    And line chart viewer should report no error
    When user hovers over the first "point" area of line chart viewer
    Then tooltip should be visible
    And tooltip should contain text "concat unique(Stereo Category)"
    And tooltip should contain text "min(Average Mass)"
    When user moves the pointer away from line chart viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown
