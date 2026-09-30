@viewers @realizes:charts.viewer.tree @realizes:viewers.filters.categorical
Feature: Tree filter panel, style settings, On Click and the Clone View and project round trips
  The Tree viewer draws a categorical hierarchy as branches. It reports no areas or readings of its
  own, so these scenarios claim that a setting reached the viewer (its property), that the viewer
  drew again (its canvas changed) and that nothing was logged; they cannot claim what a branch
  shows. The viewer is added from the ribbon's Add viewer gallery on demog (5850 rows; CONTROL true
  39) and its hierarchy is picked in the Select columns dialog of the Hierarchy property.
  Translated from the TestTrack case Charts/tree. The hierarchy is checked in the order CONTROL, SEX,
  RACE; the dialog keeps the check order, so no row is dragged. The branch clicks of charts-ui stay
  manual (see the request document).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog dataset
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Tree" card in "Add Viewer" dialog
    Then tree viewer should be visible
    And tree viewer should be bound to table "demog"
    When user clicks on grid
    And user clicks on settings icon of tree viewer
    Then "Hierarchy" property in context panel should be visible
    When user clicks on "..." button in "Hierarchy" property in context panel
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
    And tree viewer should be painted

  Scenario: The Tree follows the Filter Panel
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user takes a snapshot of tree viewer
    And user clicks on the "category true of CONTROL" area of filter panel
    Then 39 rows should pass the filter
    And the filter should pass exactly the rows where "CONTROL" is "true"
    And tree viewer should have repainted
    And no errors should have been logged
    When user takes a snapshot of tree viewer
    And user hovers over "CONTROL" filter card
    And user clicks on close of "CONTROL" filter card
    Then 5850 rows should pass the filter
    And tree viewer should have repainted
    And no errors should have been logged

  Scenario: Style, size and color settings apply without errors (github-3221, GROK-17376, GROK-17405, GROK-18087)
    Given "Style" category in context panel is expanded
    Then "Orient" property of tree viewer should be "LR"
    And "Layout" property of tree viewer should be "orthogonal"
    When user takes a snapshot of tree viewer
    And user checks "Show Counts" property in context panel
    Then "Show Counts" property of tree viewer should be "true"
    And tree viewer should have repainted
    And no errors should have been logged
    When user takes a snapshot of tree viewer
    And user enters "20" into "Font Size" property in context panel
    Then "Font Size" property of tree viewer should be "20"
    And tree viewer should have repainted
    When user takes a snapshot of tree viewer
    And user enters "0" into "Label Rotate" property in context panel
    Then "Label Rotate" property of tree viewer should be "0"
    And tree viewer should have repainted
    And no errors should have been logged
    When user takes a snapshot of tree viewer
    And user selects "TB" in "Orient" property in context panel
    Then tree viewer should have repainted
    When user takes a snapshot of tree viewer
    And user selects "RL" in "Orient" property in context panel
    Then tree viewer should have repainted
    When user takes a snapshot of tree viewer
    And user selects "LR" in "Orient" property in context panel
    Then tree viewer should have repainted
    And no errors should have been logged
    When user takes a snapshot of tree viewer
    And user selects "radial" in "Layout" property in context panel
    Then tree viewer should have repainted
    When user takes a snapshot of tree viewer
    And user selects "orthogonal" in "Layout" property in context panel
    Then tree viewer should have repainted
    And no errors should have been logged
    Given "Size" category in context panel is expanded
    When user takes a snapshot of tree viewer
    And user selects "HEIGHT" in "Size" property in context panel
    Then "Size" property of tree viewer should be "HEIGHT"
    And tree viewer should have repainted
    And no errors should have been logged
    When user selects "nulls" in "Size Aggr Type" property in context panel
    Then "Size Aggr Type" property of tree viewer should be "nulls"
    And tree viewer should be painted
    When user selects "#selected" in "Size Aggr Type" property in context panel
    Then "Size Aggr Type" property of tree viewer should be "#selected"
    And tree viewer should be painted
    When user selects "avg" in "Size Aggr Type" property in context panel
    Then "Size Aggr Type" property of tree viewer should be "avg"
    And tree viewer should be painted
    And no errors should have been logged
    Given "Color" category in context panel is expanded
    When user takes a snapshot of tree viewer
    And user selects "AGE" in "Color" property in context panel
    Then "Color" property of tree viewer should be "AGE"
    And tree viewer should have repainted
    And no errors should have been logged
    Given "Value" category in context panel is expanded
    When user unchecks "Include Nulls" property in context panel
    Then "Include Nulls" property of tree viewer should be "false"
    When user checks "Include Nulls" property in context panel
    Then "Include Nulls" property of tree viewer should be "true"
    And no errors should have been logged
    When user takes a snapshot of tree viewer
    And user double-clicks on the "cell 1 of RACE" area of grid
    And user presses Control+A
    And user types "Other" at the caret
    And user presses Enter
    Then the value of "RACE" column in row 1 should be "Other"
    And tree viewer should have repainted
    And no errors should have been logged

  Scenario: On Click sets Row Source (github-3245, GROK-18323)
    Given "Misc" category in context panel is expanded
    When user clicks on value of "On Click" property in context panel
    Then "On Click" property in context panel should offer "Select, Filter, None"
    When user selects "Filter" in "On Click" property in context panel
    Then "On Click" property of tree viewer should be "Filter"
    And "Row Source" property in context panel should contain text "All"
    When user selects "None" in "On Click" property in context panel
    Then "Row Source" property in context panel should contain text "All"
    When user selects "Select" in "On Click" property in context panel
    Then "Row Source" property in context panel should contain text "Filtered"
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 5" area of grid holding Shift
    Then 5 rows should be selected
    And rows 1 to 5 should be selected
    When user presses Escape in grid
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Settings survive Clone View and a project save and reopen (GROK-18322, GROK-18265, GROK-18324)
    Given no project named "TreeRoundTrip{time}" is on the server
    And "Style" category in context panel is expanded
    When user checks "Show Counts" property in context panel
    And user enters "20" into "Font Size" property in context panel
    And user selects "TB" in "Orient" property in context panel
    Given "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    Then properties of tree viewer should be:
      | Hierarchy   | CONTROL, SEX, RACE |
      | Show Counts | true               |
      | Font Size   | 20                 |
      | Orient      | TB                 |
      | On Click    | Filter             |
    When user picks "View > Layout > Clone View" from the top menu
    Then the open table views should be exactly "demog, demog"
    And the open tableview should have 1 tree viewer
    And tree viewer should be painted
    And properties of tree viewer should be:
      | Hierarchy   | CONTROL, SEX, RACE |
      | Show Counts | true               |
      | Font Size   | 20                 |
      | Orient      | TB                 |
      | On Click    | Filter             |
    When user closes the current view
    Then the "demog" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "TreeRoundTrip{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share TreeRoundTrip{time}" dialog should be visible
    When user clicks on CANCEL button in "Share TreeRoundTrip{time}" dialog
    Then the "Share TreeRoundTrip{time}" dialog should close
    And 1 project named "TreeRoundTrip{time}" should be on the server
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "TreeRoundTrip{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on TreeRoundTrip{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "demog" view should be current
    And tree viewer should be painted
    And properties of tree viewer should be:
      | Hierarchy   | CONTROL, SEX, RACE |
      | Show Counts | true               |
      | Font Size   | 20                 |
      | Orient      | TB                 |
      | On Click    | Filter             |
    And no errors should have been logged
