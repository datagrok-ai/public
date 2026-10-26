@viewers @realizes:charts.viewer.tree
Feature: Tree regressions of closed issues
  Closed Tree issues, each claimed by the rows a click selects or filters and by what the tree
  draws. public#3298: under On Click = Filter a Control+click adds a branch to the filter as it
  does to a selection, and after Filter the On Click = Select clicks select again (they "stopped
  working"). public#3249: rows selected in the grid are drawn on the tree (cross-viewer selection),
  and clearing the selection takes the highlight away. The viewer is added from the ribbon's Add
  viewer gallery on demog-1000 with the hierarchy CONTROL, SEX, RACE, as in tree-on-click-filter;
  false | F | Asian holds 7 rows, false | M | Asian 8.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Tree" card in "Add Viewer" dialog
    Then tree viewer should be visible
    When user clicks on grid
    And user clicks on settings icon of tree viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user toggles the "CONTROL" column in the column list of "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    And user toggles the "RACE" column in the column list of "Select columns..." dialog
    Then "3 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    And "Hierarchy" property of tree viewer should be "CONTROL, SEX, RACE"
    And tree viewer should have a "branch All | false | F | Asian" area

  Scenario: Control+click adds a branch to the filter, and Select works again after Filter (public#3298)
    Given "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    Then "On Click" property of tree viewer should be "Filter"
    When user clicks on the "branch All | false | F | Asian" area of tree viewer
    Then 7 rows should pass the filter
    When user clicks on the "branch All | false | M | Asian" area of tree viewer holding Control
    Then 15 rows should pass the filter
    And the filter should pass exactly the rows where "RACE" is "Asian"
    And no rows should be selected
    When user selects "Select" in "On Click" property in context panel
    Then "On Click" property of tree viewer should be "Select"
    When user resets the filter
    Then 1000 rows should pass the filter
    When user clicks on the "branch All | false | F | Asian" area of tree viewer
    Then 7 rows should be selected
    And no rows where "SEX" is "M" should be selected
    When user clicks on the "branch All | false | M | Asian" area of tree viewer holding Control
    Then 15 rows should be selected
    And no errors should have been logged

  Scenario: Rows selected in the grid are drawn on the tree (public#3249)
    When user takes a snapshot of tree viewer
    And user selects rows where "RACE" is "Asian"
    Then 15 rows should be selected
    And tree viewer should show more selection highlight than before
    When user takes a snapshot of tree viewer
    And user clears the row selection
    Then no rows should be selected
    And tree viewer should show less selection highlight than before
    And no errors should have been logged
