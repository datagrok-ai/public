@journey @viewers @realizes:viewers.grid
Feature: A second grid as a viewer, and column tooltips
  A grid added to a view is a viewer like any other: it shows the same table, takes its own Row
  Height and Show Column Labels without touching the view's own grid, rebinds to another table, and
  closes leaving the view's grid alone. The second one is addressed as "second grid viewer" — the
  reserved "grid" always means the table view's own. The column tooltip settings live on the header
  menu as radio items, and the item that is on says so through `aria-checked`.
  One journey on demog-1000, with spgi-100 as the second table.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And the open tableview should have 1 grid viewer

  Scenario: A second grid shows the same table and takes its own row height
    When user adds a grid viewer
    Then grid viewer should be added to the open tableview
    And the open tableview should have 2 grid viewers
    And second grid viewer should show 1000 rows
    When user sets "Row Height" property of second grid viewer to "40"
    Then the "cell 1 of USUBJID" area of second grid viewer should be at least 34 pixels tall
    And "Row Height" property of grid should not be "40"
    When user sets "Show Column Labels" property of second grid viewer to "false"
    Then second grid viewer should not have a "header AGE" area
    And grid should have a "header AGE" area
    When user sets "Show Column Labels" property of second grid viewer to "true"
    Then second grid viewer should have a "header AGE" area
    And no errors should have been logged

  Scenario: The second grid rebinds to another table and closes without taking the view's grid
    When user opens spgi dataset
    Then the table should have 100 rows
    When user switches to the "demog-1000" table view
    And user sets "Table" property of second grid viewer to "spgi-100"
    Then second grid viewer should be bound to table "spgi-100"
    And second grid viewer should show 100 rows
    And grid should show 1000 rows
    When user clicks on close icon of second grid viewer
    Then the open tableview should have 1 grid viewer
    And grid should show 1000 rows
    And grid should have a "header AGE" area
    When user closes all views
    Then no errors should have been logged

  Scenario: The column tooltip menu marks the setting that is on
    Given user opens demog-1000 dataset
    When user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    Then "Default" menu item in context menu should be visible
    And "Form" menu item in context menu should be visible
    And "Columns" menu item in context menu should be visible
    And "None" menu item in context menu should be visible
    And "Default" menu item in context menu should be selected
    And "None" menu item in context menu should not be selected
    When user closes the context menu
    And user picks "Tooltip > Current Column > None" from the context menu of the "header AGE" area of grid
    And user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    Then "None" menu item in context menu should be selected
    And "Default" menu item in context menu should not be selected
    When user clicks on "Default" menu item in context menu
    And user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    Then "Default" menu item in context menu should be selected
    And "None" menu item in context menu should not be selected
    When user closes the context menu
    Then no errors should have been logged
