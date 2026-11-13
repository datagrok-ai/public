@journey @realizes:views.projects @realizes:core#239
Feature: A dashboard saves an app view with its state and starts it only when opened for real
  core#239: an app view joins a dashboard by choice, and reopens where it was docked, on the tab, the
  sub-view and the filters it had. Usage Analysis is opened next to cars, switched to Functions >
  Execution time with the date filter set to "this month", and docked to the right of cars. The Save
  project dialog lists it left out, and its "+" includes it; the reopened dashboard brings it back to the right of cars, on Execution time, with the date filter it had. A
  Browse preview of the dashboard (a single click on it under My stuff) shows the table and does not
  start the app; pinning the preview (a double click on its tab) does.

  The Projects tab of Usage Analysis is not used: its loaders never finish on a local stand. The
  dashboard is named with the run's time (letters and digits only: the gallery search misses "-") and
  removed with its table, views and picture when the feature starts and ends.

  Background:
    Given user is logged in
    And simple mode is off
    And the user's own project "BDDLayoutApp{time}" is removed now and at feature end

  Scenario: Usage Analysis is set to a sub-view and a filter, and docked beside a table
    Given user opens cars dataset
    And user opens the "Usage Analysis" app on its "Overview" view
    When user clicks on Functions tab
    Then the "Functions" view should be current
    When user clicks on "Execution time" button
    Then "Execution time" button should be disabled
    When user types "this month" into Date text input
    And user clicks on Apply button
    Then Apply button should be disabled
    When user drags the tab of the "Functions" view to the right side of the "cars" view
    Then the "Functions" view should be docked right-of the "cars" view
    And no errors should have been logged

  Scenario: The dashboard reopens Usage Analysis where it was, on its sub-view and filter
    When user clicks on the tab of the "cars" view
    And user opens the Save project dialog from the ribbon
    Then the "Functions" view should be excluded in the Save project dialog
    When user includes the "Functions" view in the Save project dialog
    Then the "Functions" view should be included in the Save project dialog
    When user types "BDDLayoutApp{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDDLayoutApp{time}" should be on the server
    When user closes all views
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLayoutApp{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLayoutApp{time} gallery card and the project opens
    Then the "Functions" view should be docked right-of the "cars" view
    Given user switches to the "Functions" view
    Then "Execution time" button should be disabled
    And Date text input should have value "this month"
    And no errors should have been logged

  Scenario: A Browse preview of the dashboard does not start the app, and pinning it does
    When user closes all views
    Then no table should be left in the workspace
    Given the browse panel is open
    And My stuff tree node inside browse tree is expanded
    When user expands the "Dashboards" bucket of My stuff
    And user clicks on "BDDLayoutApp{time}" in the "Dashboards" bucket of My stuff
    Then "cars" view should be visible
    And "Functions" view should be absent
    When user double-clicks on cars view
    Then "Functions" view should be visible
    And no errors should have been logged
