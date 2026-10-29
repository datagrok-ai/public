@viewers @edge-cases @realizes:charts.viewer.sunburst
Feature: Sunburst edge cases: an empty hierarchy, a hierarchy column removed, On Click switched from Filter to Select
  Guards for the data shapes and the changes under an open sunburst. A hierarchy with no column left
  shows the viewer's message instead of a chart ("requires at least one categorical column",
  `sunburst-viewer.ts`). A hierarchy column removed from the table is dropped from the hierarchy
  (the viewer's column-removed handler), and the remaining rings are drawn. Switching On Click from
  Filter to Select should give back the rows a clicked segment filtered out; it does not, a known
  failure (public#4425). No rows, a renamed hierarchy column, a layout naming a missing column, and
  the filter left by a closed sunburst or by the removed column of its clicked segment are rows of
  `viewer-rules-edge-cases.feature`. Data: demog-1000 (SEX F 553 / M 447; 10 segments over SEX,
  RACE).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    And user adds a sunburst viewer with:
      | hierarchyColumnNames | SEX, RACE |
    Then the "hierarchy columns" reading of sunburst viewer should be "SEX, RACE"
    And the "segments" reading of sunburst viewer should be 10

  Scenario: A hierarchy with no column shows the viewer's message instead of a chart
    When user sets "hierarchyColumnNames" property of sunburst viewer to ""
    Then "The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories" text in sunburst viewer should be visible
    When user sets "hierarchyColumnNames" property of sunburst viewer to "SEX, RACE"
    Then "The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories" text in sunburst viewer should be absent
    And the "segments" reading of sunburst viewer should be 10
    And no errors should have been logged

  Scenario: A hierarchy column removed from the table is dropped from the hierarchy
    When user removes "SEX" column
    Then the "hierarchy columns" reading of sunburst viewer should be "RACE"
    And the "segments" reading of sunburst viewer should be 4
    And the "rows of segment Caucasian" reading of sunburst viewer should be 896
    And no errors should have been logged

  @known-failure:43 @realizes:public#4425
  Scenario: Switching On Click from Filter to Select gives back the rows its clicked segment filtered out
    When user sets "onClick" property of sunburst viewer to "Filter"
    And user clicks on the "segment F" area of sunburst viewer
    Then 553 rows should pass the filter
    When user sets "onClick" property of sunburst viewer to "Select"
    Then all rows should pass the filter
    And no errors should have been logged
