@journey @realizes:views.projects @realizes:core#239
Feature: A view docked under two views side by side reopens under those two, not under the side panels
  core#239 follow-up (a QA report): demog-1000 is docked below cars, then iris is docked to the right of
  cars, so demog-1000 runs under cars and iris both, while the Browse panel keeps the full height on the
  left. The reopened dashboard used to put demog-1000 under the Browse panel too, the whole width of the
  window: the arrangement was replayed inside-out, and iris, docked first, merged into the row the side
  panels share. It now comes back under cars and iris only.

  The dashboard is named with the run's time (letters and digits only: the gallery search misses "-") and
  removed with its tables, views and picture when the feature starts and ends.

  Background:
    Given user is logged in
    And simple mode is off
    And the user's own project "BDDLayoutSpan{time}" is removed now and at feature end

  Scenario: demog-1000 docked below cars, then iris docked to the right of cars
    Given the browse panel is open
    And user opens cars dataset
    And user opens iris dataset
    And user opens demog-1000 dataset
    When user drags the tab of the "demog-1000" view to the bottom side of the "cars" view
    Then the "demog-1000" view should be docked below the "cars" view
    When user drags the tab of the "iris" view to the right side of the "cars" view
    Then the "iris" view should be docked right-of the "cars" view
    And the "demog-1000" view should span the "cars" and "iris" views
    And no errors should have been logged

  Scenario: The reopened dashboard keeps demog-1000 under cars and iris only
    When user clicks on the tab of the "cars" view
    And user opens the Save project dialog from the ribbon
    And user types "BDDLayoutSpan{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDDLayoutSpan{time}" should be on the server
    When user closes all views
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLayoutSpan{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLayoutSpan{time} gallery card and the project opens
    Then the "iris" view should be docked right-of the "cars" view
    And the "demog-1000" view should be docked below the "cars" view
    And the "demog-1000" view should span the "cars" and "iris" views
    And no error or warning balloon should have been shown
