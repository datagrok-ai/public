@journey @viewers @realizes:viewers.stats-viewer
Feature: The Statistics and Histograms submenus, and what they turn on
  Picking an aggregation from the viewer's **Statistics** submenu adds or removes its column, and
  the viewer rebuilds `look.stats` from the grid's tagged columns — so `stats` is the record of
  what the menu did, and `header <stat>` is whether the column landed.
  The old spec right-clicked at 50 % across and 40 % down the viewer to reach the menu, read the
  item's check state out of a `fa-check|fa-dot-circle` CSS class, and proved the pick had worked by
  a five-hundred-pixel canvas diff. A checked Dart menu item now states its check through
  `aria-checked` (`widgets/menu/menu.dart:295`), so `"sum" menu item should be selected` is the
  library's generic state read with no viewer-specific binding at all, and the right-click lands on
  `row AGE`, a region the viewer reports.
  It also noted that "the submenu opens only once per page ... those steps are manual". Measured
  here: the submenu reopens, and removing `sum` again is a scenario, not a manual note.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a statistics viewer
    And user resizes statistics viewer to 900 by 400
    Then the "stats" reading of statistics viewer should be "values, nulls, unique, min, max, avg, med, stdev"
    And the "columns shown" reading of statistics viewer should be 11
    And statistics viewer should report no error

  Scenario: The Statistics submenu states which aggregations are on
    When user right-clicks on the "row AGE" area of statistics viewer
    Then the open menu should list "Statistics > sum"
    And "min" menu item should be selected
    And "max" menu item should be selected
    And "avg" menu item should be selected
    And "stdev" menu item should be selected
    And "sum" menu item should not be selected
    And "geomean" menu item should not be selected
    And "variance" menu item should not be selected
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Statistics > sum adds the column, and picking it again takes it away
    Then the "stats" reading of statistics viewer should not contain "sum"
    And statistics viewer should not have a "header sum" area
    When user picks "Statistics > sum" from the context menu of the "row AGE" area of statistics viewer
    Then the "stats" reading of statistics viewer should contain "sum"
    And statistics viewer should have a "header sum" area
    And the "sum of AGE" reading of statistics viewer should be "45677.00"
    And the "sum of SEX" reading of statistics viewer should be ""
    When user right-clicks on the "row AGE" area of statistics viewer
    Then the open menu should list "Statistics > sum"
    And "sum" menu item should be selected
    When user closes the context menu
    And user picks "Statistics > sum" from the context menu of the "row AGE" area of statistics viewer
    Then the "stats" reading of statistics viewer should not contain "sum"
    And statistics viewer should not have a "header sum" area
    And no errors should have been logged

  Scenario: An aggregation dropped from the list takes its column with it
    Then statistics viewer should have a "header med" area
    When user picks "Statistics > med" from the context menu of the "row AGE" area of statistics viewer
    Then the "stats" reading of statistics viewer should not contain "med"
    And statistics viewer should not have a "header med" area
    And statistics viewer should have a "header avg" area
    When user picks "Statistics > med" from the context menu of the "row AGE" area of statistics viewer
    Then the "stats" reading of statistics viewer should contain "med"
    And the "med of AGE" reading of statistics viewer should be "45.00"
    And no errors should have been logged

  Scenario: The Histograms submenu offers the categorical columns with fewer than ten categories
    When user right-clicks on the "row AGE" area of statistics viewer
    Then the open menu should list "Histograms > SEX"
    And the open menu should list "Histograms > RACE"
    And the open menu should list "Histograms > DIS_POP"
    And the open menu should list "Histograms > SEVERITY"
    And the open menu should not list "Histograms > USUBJID"
    And the open menu should not list "Histograms > AGE"
    And the open menu should not list "Histograms > STARTED"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Histograms > SEX adds a histogram column, and picking it again removes it
    Then the "histogram columns" reading of statistics viewer should be 0
    When user picks "Histograms > SEX" from the context menu of the "row AGE" area of statistics viewer
    Then the "histogram columns" reading of statistics viewer should be 1
    And the "stats" reading of statistics viewer should contain "avg"
    And the "stats" reading of statistics viewer should contain "med"
    And the "stats" reading of statistics viewer should contain "stdev"
    And the "stats" reading of statistics viewer should not contain "sum"
    When user picks "Histograms > SEX" from the context menu of the "row AGE" area of statistics viewer
    Then the "histogram columns" reading of statistics viewer should be 0
    And no errors should have been logged
