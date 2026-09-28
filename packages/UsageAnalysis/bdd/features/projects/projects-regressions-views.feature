@serial @realizes:views.projects
Feature: Projects regressions: the views a reopened project shows
  Regression guards for fixed Jira bugs about which views a project reopens with, and how they are
  arranged: GROK-13303 (the tab order), GROK-18454 (the active view), GROK-18607 (the view selector
  with tabs turned off) and GROK-17765 (reopening after each view was closed by its tab). Every
  project is saved through the ribbon's Save dialog and reopened from the Dashboards gallery. The
  demo tables are opened from Browse > Files > Demo; demog, cars and iris stand for the four tables
  of GROK-13303 and the SPGI project of GROK-18607. As in GROK-13303, the tables are renamed
  (Table > Rename... of each view, to TabOne, TabTwo and TabThree) before the tabs are dragged and
  the project is saved: a view whose name no longer matches the saved layout is the one the reopen
  puts out of order.

  Each scenario fails on the behaviour its ticket describes: an order other than the dragged one
  (the opening order, its reverse or alphabetical), the first or the last tab active instead of the
  middle one, no view selector, the "already open" warning with one view.

  The same two claims are made once more after a second save of the reopened project ("Save
  original project"), which must replace the first save: one project of that name stays on the
  server, and it reopens with the new layout.

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after a save lets
  that one message through; every check after a reopen is strict.

  Every project is named with the run's time and removed with its tables and views when each
  scenario starts and when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-13303
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
    When user clicks on browse tab
    Then the browse tree should be visible
    Then the table view tabs should read "iris, cars, demog"
    When user clicks on iris view
    Then the "iris" view should be current
    When user picks "Table > Rename..." from the context menu of iris view
    Then "Rename table" dialog should be visible
    When user enters "TabThree" into "New name:" text input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And the "TabThree" view should be current
    When user clicks on cars view
    Then the "cars" view should be current
    When user picks "Table > Rename..." from the context menu of cars view
    Then "Rename table" dialog should be visible
    When user enters "TabTwo" into "New name:" text input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And the "TabTwo" view should be current
    When user clicks on demog view
    Then the "demog" view should be current
    When user picks "Table > Rename..." from the context menu of demog view
    Then "Rename table" dialog should be visible
    When user enters "TabOne" into "New name:" text input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And the "TabOne" view should be current
    Then the table view tabs should read "TabThree, TabTwo, TabOne"
    And the open tables should be exactly "TabOne, TabTwo, TabThree"
    When user drags TabTwo view to TabOne view
    Then the table view tabs should read "TabThree, TabOne, TabTwo"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegTabs{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegTabs{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegTabs{time}" dialog
    Then the "Share BDDRegTabs{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegTabs{time}" into gallery search
    And user double-clicks on BDDRegTabs{time} gallery card
    Then the table view tabs should read "TabThree, TabOne, TabTwo"
    And no errors should have been logged

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
    Then the table view tabs should read "iris, cars, demog"
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
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegActive{time}" into gallery search
    And user double-clicks on BDDRegActive{time} gallery card
    Then the table view tabs should read "iris, cars, demog"
    And the "cars" view should be current
    And cars view should be selected
    And no errors should have been logged

  @realizes:GROK-13303 @realizes:GROK-18454
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
    When user clicks on browse tab
    Then the browse tree should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegResave{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegResave{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegResave{time}" dialog
    Then the "Share BDDRegResave{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegResave{time}" into gallery search
    And user double-clicks on BDDRegResave{time} gallery card
    Then the table view tabs should read "iris, cars, demog"
    And the "iris" view should be current
    When user drags cars view to demog view
    And user clicks on demog view
    Then the table view tabs should read "iris, demog, cars"
    And the "demog" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegResave{time}" uploaded' should have been shown
    And 1 project named "BDDRegResave{time}" should be on the server
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegResave{time}" into gallery search
    And user double-clicks on BDDRegResave{time} gallery card
    Then the table view tabs should read "iris, demog, cars"
    And the "demog" view should be current
    And no errors should have been logged

  @realizes:GROK-18607
  Scenario: With the view tabs turned off, a reopened project shows the view selector
    Given no project named "BDDRegNoTabs{time}" is on the server
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
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegNoTabs{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegNoTabs{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegNoTabs{time}" dialog
    Then the "Share BDDRegNoTabs{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Tabs button inside status bar
    Then the view tabs should be hidden
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegNoTabs{time}" into gallery search
    And user double-clicks on BDDRegNoTabs{time} gallery card
    Then table "demog" should be open
    And table "cars" should be open
    And the view selector should list "demog, cars"
    And no errors should have been logged

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
    When user clicks on browse tab
    Then the browse tree should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegReopen{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegReopen{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegReopen{time}" dialog
    Then the "Share BDDRegReopen{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegReopen{time}" into gallery search
    And user double-clicks on BDDRegReopen{time} gallery card
    Then the table view tabs should read "cars, demog"
    When user closes the "demog" view by the cross on its tab
    And user closes the "cars" view by the cross on its tab
    Then the table view tabs should read ""
    When user double-clicks on BDDRegReopen{time} gallery card
    Then the table view tabs should read "cars, demog"
    And no error or warning balloon should have been shown
    And no errors should have been logged
