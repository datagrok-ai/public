@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — a configured cross tab that survives a round trip
  A pivot configured away from its defaults — med instead of avg, a title, a colour-coded value
  column — saved with the view's layout on the server and restored, and then saved as a project,
  closed and reopened: the tag rows, the aggregation it recomputes and the inner grid's look all come
  back (github-2535). One journey on demog-1000 grouped by DIS_POP and pivoted on SEVERITY, where
  med(AGE) is a different number from avg(AGE) in every cell, so a restored viewer that fell back to
  the default aggregation would be caught.
  Not translated: the ribbon's Save entry point of the old server spec — generic application chrome,
  and the project round-trip below already proves the channel.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    And user picks "Aggregation > med" from the context menu of the "aggregate chip avg(AGE)" area of pivot table viewer
    And user closes the context menu
    Then the "aggregations" reading of pivot table viewer should be "med(AGE)"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the aggregated values of pivot table viewer should match "med(AGE)" grouped by "DIS_POP" pivoted on "SEVERITY"

  Scenario: A layout saved on the server brings the non-default aggregation back
    When user saves the layout of the current table view to the server
    And user clicks on close icon of pivot table viewer
    Then pivot table viewer should be absent
    When user loads the saved layout
    Then pivot table viewer should be visible
    And the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "pivot" reading of pivot table viewer should be "SEVERITY"
    And the "aggregate" reading of pivot table viewer should be "med(AGE)"
    And the "aggregations" reading of pivot table viewer should be "med(AGE)"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the aggregated values of pivot table viewer should match "med(AGE)" grouped by "DIS_POP" pivoted on "SEVERITY"
    And no errors should have been logged

  Scenario: The title and the inner grid's colour coding travel with the layout
    When user sets "Title" property of pivot table viewer to "Cross tab"
    And user picks "Grid > Color Coding > Linear" from the context menu of the "grid header None med(AGE)" area of pivot table viewer
    Then the "color coding of None med(AGE)" reading of pivot table viewer should be "Linear"
    And title of pivot table viewer should have text "Cross tab"
    When user saves the layout of the current table view to the server
    And user clicks on close icon of pivot table viewer
    Then pivot table viewer should be absent
    When user loads the saved layout
    Then pivot table viewer should be visible
    And title of pivot table viewer should have text "Cross tab"
    And the "color coding of None med(AGE)" reading of pivot table viewer should be "Linear"
    And the "color of grid cell 2 of None med(AGE)" reading of pivot table viewer should not be "#ffffff"
    And the "aggregations" reading of pivot table viewer should be "med(AGE)"
    And no errors should have been logged

  Scenario: A project saved, closed and reopened brings the whole configuration back
    When user saves the current view as project "bdd pivot table round trip"
    And user closes all views
    And user opens the "bdd pivot table round trip" project
    Then pivot table viewer should be visible
    And the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "pivot" reading of pivot table viewer should be "SEVERITY"
    And the "aggregations" reading of pivot table viewer should be "med(AGE)"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "aggregated columns" reading of pivot table viewer should be 6
    And title of pivot table viewer should have text "Cross tab"
    And the aggregated values of pivot table viewer should match "med(AGE)" grouped by "DIS_POP" pivoted on "SEVERITY"
    And no errors should have been logged
