@viewers @realizes:charts.viewer.tree
Feature: Tree clicks: On Click = Filter and Select, Shift and Control, and the grid's selection drawn on the tree
  With On Click set to Filter, a click on a branch of the Tree filters the table to the rows under it,
  a Shift+click or a Control+click adds another branch's rows, and a plain click starts over from the
  new branch; switched back to Select, the clicks select again (public#3298). Rows selected in the
  grid are drawn on the tree, and clearing the selection takes the highlight away (public#3249).
  tree.feature claims that the property switches Row Source; here the clicks themselves are claimed.
  The viewer is added from the ribbon's Add viewer gallery on demog-1000 and its hierarchy is
  CONTROL, SEX, RACE, picked in the Select columns dialog. Counts on demog-1000: false | F | Asian 7,
  false | M | Asian 8 (every Asian row), true | F | Caucasian 5. A branch is the viewer's
  `branch <path>` area, the line into the node, where a click lands on the node it leads into.
  A first click made with Shift filters to the branch too (public#4393, fixed).
  Translated from the TestTrack case Charts/tree (On Click) and the github-3245 spec
  tree-rowsource-onclick-state-bug.

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

  @realizes:public#3249 @realizes:public#3298
  Scenario: The grid's selection is drawn on the tree; a click filters to the branch, Shift or Control adds a branch, a plain click starts over, and Select selects again (public#3249, public#3298)
    When user takes a snapshot of tree viewer
    And user selects rows where "RACE" is "Asian"
    Then 15 rows should be selected
    And tree viewer should show more selection highlight than before
    When user takes a snapshot of tree viewer
    And user clears the row selection
    Then no rows should be selected
    And tree viewer should show less selection highlight than before
    Given "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    Then "On Click" property of tree viewer should be "Filter"
    And 1000 rows should pass the filter
    When user clicks on the "branch All | false | F | Asian" area of tree viewer
    Then 7 rows should pass the filter
    And no rows where "SEX" is "M" should pass the filter
    And no rows where "RACE" is "Caucasian" should pass the filter
    And no rows should be selected
    When user clicks on the "branch All | false | M | Asian" area of tree viewer holding Shift
    Then 15 rows should pass the filter
    And the filter should pass exactly the rows where "RACE" is "Asian"
    When user clicks on the "branch All | true | F | Caucasian" area of tree viewer
    Then 5 rows should pass the filter
    And no rows where "CONTROL" is "false" should pass the filter
    And no rows where "RACE" is "Asian" should pass the filter
    And no rows should be selected
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

  Scenario: A first click made with Shift filters to the branch as a plain click does
    Given "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    Then "On Click" property of tree viewer should be "Filter"
    And 1000 rows should pass the filter
    When user clicks on the "branch All | false | F | Asian" area of tree viewer holding Shift
    Then 7 rows should pass the filter
    And no errors should have been logged
