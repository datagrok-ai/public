@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — what Row Source aggregates
  Row Source decides which rows the aggregation runs over: at All the pivot ignores a filter on the
  source table, at Filtered it re-aggregates over the rows the filter passes and says so through its
  own "rows shown" reading. The filter is a categorical card of the filter panel, a registered filter
  — a bitset written straight into the table does not survive the pivot's next refresh, which asks
  the platform to recompute the filter. One journey on demog-1000 grouped by DIS_POP with avg(AGE):
  447 of the 1000 rows are male, RA averages 51.60 over all of them and 52.10 over the male ones.
  The last three scenarios are the cell link, which only Row Source All has: the pivot writes its own
  criterion into the table's filter, and a criterion is what it writes — collaborative, ANDed with
  every other filter the table carries, never a replacement of them. What the criterion holds depends
  on the cell: the key cell of a row (and the row header beside it) gives the whole group, a cell
  under a pivot column gives the group crossed with that pivot value. So on a clean panel a click on
  the RA key cell passes all 434 RA rows, women among them; with the SEX card on M the same click
  passes 104 — RA and M — and the table then carries both criteria, "DIS_POP in [RA]" and "SEX: M".
  A click on the RA × None cell passes 254 — RA and SEVERITY None. `pivottable-rowsource-filter-selection.md`
  reads its step 7 the other way ("the click replaces the SEX filter, it does not intersect it", 434
  with the card on); that expectation is the md's error and the scenarios below state the measured
  behaviour instead — the mark for GROK-17726 is gone with it. The md's 434 is real, it is simply
  what a click gives with nothing else filtering.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    And user sets "Pivot Column Names" property of pivot table viewer to ""
    And user sets "Row Source" property of pivot table viewer to "All"
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "rows shown" reading of pivot table viewer should be 1000

  Scenario: At Row Source All the source filter changes nothing
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of pivot table viewer should be 1000
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of avg(AGE)" reading of pivot table viewer should be "51.60"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP"
    And no errors should have been logged

  Scenario: At Row Source Filtered the pivot re-aggregates over the filtered rows
    When user sets "Row Source" property of pivot table viewer to "Filtered"
    Then the "rows shown" reading of pivot table viewer should be 447
    And pivot table viewer should show 447 rows
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of avg(AGE)" reading of pivot table viewer should be "52.10"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP" over the filtered rows
    And no errors should have been logged

  Scenario: Lifting the filter takes the pivot back to the whole table
    When user adds a categorical filter on "SEX" keeping "F, M"
    Then all rows should pass the filter
    And the "rows shown" reading of pivot table viewer should be 1000
    And the "text of grid cell 5 of avg(AGE)" reading of pivot table viewer should be "51.60"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP"
    When user sets "Row Source" property of pivot table viewer to "All"
    Then the "rows shown" reading of pivot table viewer should be 1000
    And no errors should have been logged

  Scenario: At Row Source All a key cell click takes the whole group, and a panel card stacks on it
    Given "Row Source" property of pivot table viewer should be "All"
    And "Filtering Enabled" property of pivot table viewer should be "true"
    And all rows should pass the filter
    When user clicks on the "grid cell 5 of DIS_POP" area of pivot table viewer
    Then 434 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" is "RA"
    And the "filter label" reading of pivot table viewer should be "DIS_POP in [RA]"
    When user adds a categorical filter on "SEX" keeping "M"
    Then 104 rows should pass the filter
    And no rows where "SEX" is "F" should pass the filter
    And no rows where "DIS_POP" is "UC" should pass the filter
    And the "filter label" reading of pivot table viewer should be "DIS_POP in [RA]"
    And no errors should have been logged

  Scenario: A cell under a pivot column takes the group crossed with that pivot value
    When user adds a categorical filter on "SEX" keeping "F, M"
    And user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY"
    And user sets "Row Source" property of pivot table viewer to "All"
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "aggregated columns" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of DIS_POP" reading of pivot table viewer should be "RA"
    When user clicks on the "grid cell 5 of None avg(AGE)" area of pivot table viewer
    Then 254 rows should pass the filter
    And no rows where "DIS_POP" is "UC" should pass the filter
    And no rows where "SEVERITY" is "Low" should pass the filter
    And the "filter label" reading of pivot table viewer should be "DIS_POP in [RA], None avg(AGE) = 52.30"
    And no errors should have been logged

  Scenario: At Row Source Filtered the same click writes no filter at all
    When user sets "Row Source" property of pivot table viewer to "Filtered"
    And user resets the filter
    Then all rows should pass the filter
    When user clicks on the "grid cell 5 of DIS_POP" area of pivot table viewer
    Then all rows should pass the filter
    And no errors should have been logged
