@journey @viewers @realizes:charts.viewer.timelines
Feature: Timelines regressions of closed issues
  Closed Timelines issues, each claimed by what the viewer draws. public#617: Color set and then
  cleared (none) threw an exception; it now draws every one of the 143 intervals of ae.csv with no
  legend and logs nothing. public#1973: a Split By column whose categories are very long showed no
  visualization; a copy of USUBJID with a long tail on every value splits into 71 lanes, one per
  subject as USUBJID does, with all 143 intervals drawn. The calculated column stays in the
  feature's own copy of ae, which the shell reset closes. ae.csv is opened from Browse and the
  viewer added from the ribbon's Add viewer gallery, as in timelines.feature.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Charts tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---Charts---ae.csv tree node inside browse tree
    Then the "ae" view should be current
    And the table should have 143 rows
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Timelines" card in "Add Viewer" dialog
    Then timelines viewer should be visible
    And timelines viewer should be painted
    And the "intervals" reading of timelines viewer should be 143
    When user clicks on grid
    And user clicks on settings icon of timelines viewer

  @realizes:public#617
  Scenario: Clearing Color keeps every interval and throws nothing (public#617)
    Given "Color" category in context panel is expanded
    When user selects "AESOC" in "Color" property in context panel
    Then the legend of timelines viewer should list 15 items
    When user takes a snapshot of timelines viewer
    And user sets "colorColumnName" property of timelines viewer to ""
    Then timelines viewer should have repainted
    And legend of timelines viewer should be hidden
    And the "intervals" reading of timelines viewer should be 143
    And timelines viewer should report no error
    And no errors should have been logged

  @realizes:public#1973
  Scenario: Long Split By categories still draw their lanes (public#1973)
    Given "Data" category in context panel is expanded
    Then "Split By" property of timelines viewer should be "USUBJID"
    And the "lanes" reading of timelines viewer should be 71
    When user adds a calculated column "LongSubj" with formula "${USUBJID} + ' is a deliberately very long subject label that runs past the axis'"
    And user takes a snapshot of timelines viewer
    And user selects "LongSubj" in "Split By" property in context panel
    Then "Split By" property of timelines viewer should be "LongSubj"
    And timelines viewer should have repainted
    And the "lanes" reading of timelines viewer should be 71
    And the "intervals" reading of timelines viewer should be 143
    And timelines viewer should report no error
    And timelines viewer should be painted
    When user selects "USUBJID" in "Split By" property in context panel
    Then "Split By" property of timelines viewer should be "USUBJID"
    And no errors should have been logged
