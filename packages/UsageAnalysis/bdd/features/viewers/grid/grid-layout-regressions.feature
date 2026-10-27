@viewers @realizes:viewers.grid
Feature: Grid regressions — what the grid keeps after a layout is applied
  Closed grid issues that showed only after a saved layout was applied to the view: the pinned columns
  are drawn side by side, not on top of each other, when a viewer is stacked over the grid as a tab
  (public#1625, reported on SPGI with Id and Structure; demog-1000's USUBJID and AGE need no package
  renderer), and a click on a row number still makes that row current and a Shift-click on another
  selects the rows between them (public#679). Each scenario saves its own layout to the server, and the
  layouts go at feature end.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  # public#1625
  @realizes:public#1625
  Scenario: Pinned columns stay side by side after a layout with a viewer stacked over the grid (public#1625)
    When user picks "Pin > Pin Column" from the context menu of the "header USUBJID" area of grid
    And user picks "Pin > Pin Column" from the context menu of the "header AGE" area of grid
    Then the grid should pin the columns "USUBJID, AGE"
    And the "header AGE" area of grid should lie to the right of the "header USUBJID" area
    Given user adds a scatter plot viewer
    When user stacks scatter plot viewer onto grid as a tab
    Then the viewers "demog-1000, Scatter plot" should be the tabs of one panel
    When user saves the layout of the current table view to the server
    And user loads the saved layout
    And user switches the tabbed panel to the "demog-1000" tab
    Then the grid should pin the columns "USUBJID, AGE"
    And the "header AGE" area of grid should lie to the right of the "header USUBJID" area
    And no errors should have been logged

  # public#679
  @realizes:public#679
  Scenario: A row-number click still makes the row current after a layout is applied (public#679)
    Given user adds a scatter plot viewer
    When user saves the layout of the current table view to the server
    And user loads the saved layout
    Then scatter plot viewer should be visible
    When user clicks on the "row header 5" area of grid
    Then row 5 should be current
    When user clicks on the "row header 7" area of grid holding Shift
    Then rows 5 to 7 should be selected
    And 3 rows should be selected
    And no errors should have been logged
