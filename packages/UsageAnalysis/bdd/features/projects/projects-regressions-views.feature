@serial @realizes:views.projects
Feature: Projects regressions: the views a reopened project brings back
  Regression guard for GROK-18454 (the active view): demog, cars and iris are opened from Browse >
  Files > Demo, the middle tab (cars) is made current, the project is saved through the ribbon's Save
  dialog and reopened from the Dashboards gallery. The scenario fails on the behaviour the ticket
  describes: the first or the last tab active instead of the middle one.

  public#2060 (the tab order a reopened project keeps, also after a second save, with the active
  view of GROK-18454) and GROK-17765 (reopening after each view was closed by the cross on its tab)
  are read in the tab strip of the full shell (simple mode off): the order of the views' tabs, and
  each tab gone after its cross. A tab dropped on another lands right after it; new tabs open on
  the left. Parked in the request document: GROK-18607 (the view selector with the view tabs turned
  off), which needs a reading of the view selector.

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

  @realizes:public#2060
  Scenario: A reopened project keeps the order its view tabs were dragged into
    Given no project named "BDDRegTabs{time}" is on the server
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
    And the view tabs should be in the order "iris, cars, demog"
    When user drags "cars" view tab to "demog" view tab
    Then the view tabs should be in the order "iris, demog, cars"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegTabs{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegTabs{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegTabs{time}" dialog
    Then the "Share BDDRegTabs{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Then no table should be left in the workspace
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegTabs{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegTabs{time} gallery card
    Then the view tabs should be in the order "iris, demog, cars"
    And table "cars" should be open
    And no error or warning balloon should have been shown

  @realizes:public#2060 @realizes:GROK-18454
  Scenario: A reopened project saved again comes back with its new tab order and active view
    Given no project named "BDDRegResave{time}" is on the server
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
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegResave{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegResave{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegResave{time}" dialog
    Then the "Share BDDRegResave{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Then no table should be left in the workspace
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegResave{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegResave{time} gallery card
    Then the view tabs should be in the order "iris, cars, demog"
    And the "iris" view should be current
    When user drags "cars" view tab to "demog" view tab
    And user clicks on demog view
    Then the view tabs should be in the order "iris, demog, cars"
    And the "demog" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegResave{time}" uploaded' should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Then no table should be left in the workspace
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegResave{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegResave{time} gallery card
    Then the view tabs should be in the order "iris, demog, cars"
    And the "demog" view should be current
    And no error or warning balloon should have been shown

  @realizes:GROK-17765
  Scenario: A project whose views were closed one by one opens again with all of them
    Given no project named "BDDRegReopen{time}" is on the server
    And simple mode is off
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegReopen{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegReopen{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegReopen{time}" dialog
    Then the "Share BDDRegReopen{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Then no table should be left in the workspace
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegReopen{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegReopen{time} gallery card
    Then the view tabs should be in the order "cars, demog"
    When user closes "demog" view tab
    Then "demog" view tab should be absent
    When user closes "cars" view tab
    Then "demog" view tab should be absent
    And "cars" view tab should be absent
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegReopen{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegReopen{time} gallery card
    Then the view tabs should be in the order "cars, demog"
    And no error or warning balloon should have been shown
