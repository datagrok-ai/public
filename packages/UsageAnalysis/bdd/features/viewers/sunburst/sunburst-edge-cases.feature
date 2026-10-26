@viewers @edge-cases @realizes:charts.viewer.sunburst
Feature: Sunburst edge cases: no rows, an empty hierarchy, a hierarchy column removed, a layout with a missing column
  Guards for the data shapes and the changes under an open sunburst. A filter that removes every row
  leaves nothing to draw and no error, and the rings come back with the rows. A hierarchy with no
  column left shows the viewer's message instead of a chart ("requires at least one categorical
  column", `sunburst-viewer.ts`). A hierarchy column removed from the table is dropped from the
  hierarchy (the viewer's column-removed handler), and the remaining rings are drawn. In Filter mode
  the filter a clicked segment set is cleared by a double-click on empty space (help sunburst.md,
  "Clear the sunburst's filter"), also after the clicked segment's column was removed. Whether the
  removal itself should release that filter is not defined: see MISSING.md. A layout naming a
  hierarchy column the table no longer has is applied without an error. Closing the sunburst should
  give back the rows its clicked segment filtered out (a filter that outlives its viewer is a defect,
  core#233); it does not, a known failure with no ticket yet (`@ticket-pending`). Data: demog-1000
  (SEX F 553 / M 447; 10 segments over SEX, RACE).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    And user adds a sunburst viewer with:
      | hierarchyColumnNames | SEX, RACE |
    Then the "hierarchy columns" reading of sunburst viewer should be "SEX, RACE"
    And the "segments" reading of sunburst viewer should be 10

  Scenario: A filter that removes every row leaves the sunburst empty without an error
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And the "rows shown" reading of sunburst viewer should be 0
    And sunburst viewer should not have a "segment F" area
    When user resets the filter
    Then the "rows shown" reading of sunburst viewer should be 1000
    And the "segments" reading of sunburst viewer should be 10
    And no errors should have been logged

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

  Scenario: A double-click on empty space clears the filter after the clicked segment's column is removed
    When user sets "onClick" property of sunburst viewer to "Filter"
    And user clicks on the "segment F" area of sunburst viewer
    Then 553 rows should pass the filter
    When user removes "SEX" column
    Then the "hierarchy columns" reading of sunburst viewer should be "RACE"
    When user double-clicks on empty plot space of sunburst viewer
    Then 1000 rows should pass the filter
    And no errors should have been logged

  Scenario: A layout naming a hierarchy column the table no longer has is applied without an error
    When user saves the layout of the current table view
    And user removes "SEX" column
    And user loads the saved layout
    Then sunburst viewer should be visible
    And the "hierarchy columns" reading of sunburst viewer should be "RACE"
    And no errors should have been logged

  @known-failure:77 @ticket-pending
  Scenario: Closing the sunburst gives back the rows its clicked segment filtered out
    Closing the sunburst leaves the table filtered to the clicked segment, with no viewer left to
    clear it. No ticket yet. Remove the tag when closing the sunburst gives every row back.
    When user sets "onClick" property of sunburst viewer to "Filter"
    And user clicks on the "segment F" area of sunburst viewer
    Then 553 rows should pass the filter
    When user clicks on close icon of sunburst viewer
    Then the open tableview should have 0 sunburst viewers
    And all rows should pass the filter
    And no errors should have been logged
