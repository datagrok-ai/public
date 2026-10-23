@viewers @realizes:charts.viewer.tree
Feature: Tree with On Click = Filter filters the table to the clicked branch
  With On Click set to Filter, a click on a branch of the Tree filters the table to the rows under it,
  a Shift+click adds another branch's rows, and a plain click starts over from the new branch.
  tree.feature claims that the property switches Row Source; here the clicks themselves are claimed.
  The viewer is added from the ribbon's Add viewer gallery on demog and its hierarchy is CONTROL,
  SEX, RACE, picked in the Select columns dialog. Counts on demog (5850 rows): false | F | Asian 37,
  false | M | Asian 35, true | F | Black 2. A branch is the viewer's `branch <path>` area, the line
  into the node, where a click lands on the node it leads into.
  Translated from the TestTrack case Charts/tree (On Click) and the github-3245 spec
  tree-rowsource-onclick-state-bug.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog dataset
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

  Scenario: A click filters to the branch, Shift+click adds a branch, a plain click starts over
    Given "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    Then "On Click" property of tree viewer should be "Filter"
    And 5850 rows should pass the filter
    When user clicks on the "branch All | false | F | Asian" area of tree viewer
    Then 37 rows should pass the filter
    And no rows where "SEX" is "M" should pass the filter
    And no rows where "RACE" is "Caucasian" should pass the filter
    And no rows should be selected
    When user clicks on the "branch All | false | M | Asian" area of tree viewer holding Shift
    Then 72 rows should pass the filter
    And the filter should pass exactly the rows where "RACE" is "Asian"
    When user clicks on the "branch All | true | F | Black" area of tree viewer
    Then 2 rows should pass the filter
    And no rows where "CONTROL" is "false" should pass the filter
    And no rows where "RACE" is "Asian" should pass the filter
    And no rows should be selected
    And no errors should have been logged
