@realizes:views.projects
Feature: Projects regressions: the views a reopened project brings back
  Regression guards for GROK-18454 (the active view) and public#2060 (the tab order) on one
  project: demog, cars and iris are opened from Browse > Files > Demo, saved through the ribbon's
  Save dialog and reopened from the Dashboards gallery; cars is then dragged after demog, the middle
  tab (demog) made current, and the project saved again and reopened. GROK-18454 is the middle tab
  coming back current, not the first or the last; public#2060 the order the drag left. GROK-17765
  (reopening after each view was closed by the cross on its tab) is the second scenario. Both are
  read in the tab strip of the full shell (simple mode off): the order of the views' tabs, and each
  tab gone after its cross. A tab dropped on another lands right after it; new tabs open on the left.
  Parked in the request document: GROK-18607 (the view selector with the view tabs turned off), which
  needs a reading of the view selector.

  The error floor is not claimed: the Save dialog's preview logs "Unable to find element in cloned
  iframe" for these views on every run (GROK-18606, won't fix), and the check that lets only that
  message through is requested in the request document.

  The project is named with the run's time and removed with its tables and views when the scenario
  starts and when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open

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
    And user double-clicks on BDDRegResave{time} gallery card and the project opens
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
    And user double-clicks on BDDRegResave{time} gallery card and the project opens
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
    And user double-clicks on BDDRegReopen{time} gallery card and the project opens
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
    And user double-clicks on BDDRegReopen{time} gallery card and the project opens
    Then the view tabs should be in the order "cars, demog"
    And no error or warning balloon should have been shown
