@viewers @realizes:viewers.pivot-table
Feature: Pivot table — the Table property switched to another table and back
  "Scenario 8 Step 7" of `pivot-table.md` (github-3414, GROK-14995): the Table property switched to
  another table and back leaves the Data row with that one table, not a second copy of it, and the
  pivot aggregates the table it is bound to each time. Claimed on the Data row's chips (`data`, the
  captions it draws, comma-joined) and on the cross tab itself, where the old spec counted the
  row's DOM entries.
  Not translated: the tag-editor header that GROK-14995 duplicated — nothing reports how many header
  rows the panel holds, see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens spgi dataset
    And user switches to the "demog-1000" table view
    And user adds a pivot table viewer
    Then the "data" reading of pivot table viewer should be "demog-1000"
    And the "aggregated rows" reading of pivot table viewer should be 6

  Scenario: Switched to spgi-100 and back, the Data row names one table each time and the cross tab follows
    When user sets "table" property of pivot table viewer to "spgi-100"
    Then pivot table viewer should be bound to table "spgi-100"
    And the "data" reading of pivot table viewer should be "spgi-100"
    And the "rows shown" reading of pivot table viewer should be 100
    And the "error" reading of pivot table viewer should be ""
    When user sets "table" property of pivot table viewer to "demog-1000"
    Then pivot table viewer should be bound to table "demog-1000"
    And the "data" reading of pivot table viewer should be "demog-1000"
    And the "rows shown" reading of pivot table viewer should be 1000
    And the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And pivot table viewer should have a "data chip demog-1000" area
    And pivot table viewer should not have a "data chip spgi-100" area
    And no errors should have been logged
