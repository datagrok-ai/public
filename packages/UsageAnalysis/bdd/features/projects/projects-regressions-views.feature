@serial @realizes:views.projects
Feature: Projects regressions: the view a reopened project makes current
  Regression guard for GROK-18454 (the active view): demog, cars and iris are opened from Browse >
  Files > Demo, the middle tab (cars) is made current, the project is saved through the ribbon's Save
  dialog and reopened from the Dashboards gallery. The scenario fails on the behaviour the ticket
  describes: the first or the last tab active instead of the middle one.

  Parked in the request document until the library has the phrases: GROK-13303 (the tab
  order a reopened project keeps, also after a second save), GROK-18607 (the view selector with the
  view tabs turned off) and GROK-17765 (reopening after each view was closed by the cross on its
  tab), which read the order of the view tabs, the view selector and a view tab's cross.

  The error floor is not claimed: the Save dialog's preview logs "Unable to find element in cloned
  iframe" for these views on every run (GROK-18606, won't fix), and the check that lets only that
  message through is requested in the request document.

  The project is named with the run's time and removed with its tables and views when the scenario
  starts and when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-18454
  Scenario: A reopened project makes current the view that was current when it was saved
    Given no project named "BDDRegActive{time}" is on the server
    And simple mode is off
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user double-clicks Files---Demo---iris.csv tree node inside browse tree
    Then the "iris" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user clicks on cars view
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegActive{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegActive{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegActive{time}" dialog
    Then the "Share BDDRegActive{time}" dialog should close
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegActive{time}" into gallery search
    And user double-clicks on BDDRegActive{time} gallery card
    Then demog view should be visible
    And iris view should be visible
    And the "cars" view should be current
    And cars view should be selected
