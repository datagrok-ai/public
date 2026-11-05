@viewers @realizes:charts.viewer.timelines
Feature: Timelines legend, legend visibility, legend clicks and Reset View
  The Timelines viewer draws one lane per subject (Split By) with an interval per event, colored by
  Color, and reports what it drew: `lanes`, `intervals` (after the legend and the zoom) and a `view`
  area over its plot. Its legend is the standard DOM legend, so its items are counted and clicked by
  name. ae.csv (143 adverse events, 71 subjects; AESOC has 15 values) is opened from Browse and the
  viewer is added from the ribbon's Add viewer gallery. Color offers the string columns only.
  Translated from the TestTrack case Charts/timelines.

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
    And timelines viewer should be bound to table "ae"
    And timelines viewer should be painted
    When user clicks on grid
    And user clicks on settings icon of timelines viewer
    Given "Color" category in context panel is expanded
    Then "Color" property in context panel should be visible

  Scenario: The viewer draws, and the legend follows Legend Visibility (GROK-20800)
    Then no errors should have been logged
    And the "intervals" reading of timelines viewer should be 143
    And "Color" property in context panel should offer the columns "STUDYID, DOMAIN, USUBJID, AESPID, AETERM, AELLT, AELLTCD, AEDECOD, AEPTCD, AEHLT, AEHLTCD, AEHLGT, AEHLGTCD, AEBODSYS, AEBDSYCD, AESOC, AESOCCD, AESEV, AEACN, AEREL, AEOUT, AESTDTC"
    And "Color" property in context panel should not offer the column "AESTDY"
    And "Color" property in context panel should not offer the column "AESEQ"
    When user selects "AESOC" in "Color" property in context panel
    Then "Color" property of timelines viewer should be "AESOC"
    And the legend of timelines viewer should list 15 items
    And legend of timelines viewer should be visible
    Given "Legend" category in context panel is expanded
    When user selects "Always" in "Legend Visibility" property in context panel
    Then legend of timelines viewer should be visible
    And the legend of timelines viewer should list 15 items
    When user selects "Never" in "Legend Visibility" property in context panel
    Then legend of timelines viewer should be hidden
    When user selects "Auto" in "Legend Visibility" property in context panel
    Then legend of timelines viewer should be visible
    And the legend of timelines viewer should list 15 items
    And no errors should have been logged

  Scenario: Legend clicks narrow the viewer without blanking it (GROK-19033, GROK-19535, GROK-18608)
    When user selects "AESOC" in "Color" property in context panel
    Then the legend of timelines viewer should list 15 items
    And the "intervals" reading of timelines viewer should be 143
    When user clicks on "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" item in the legend of timelines viewer
    Then the "intervals" reading of timelines viewer should be 35
    And 143 rows should pass the filter
    And no errors should have been logged
    When user clicks on "NERVOUS SYSTEM DISORDERS" item in the legend of timelines viewer holding Control
    Then the "intervals" reading of timelines viewer should be 68
    And no errors should have been logged
    When user clicks on "NERVOUS SYSTEM DISORDERS" item in the legend of timelines viewer holding Control
    Then the "intervals" reading of timelines viewer should be 35
    When user clicks on "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" item in the legend of timelines viewer
    Then the "intervals" reading of timelines viewer should be 143
    And "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" legend item in legend of timelines viewer should not be selected
    And no errors should have been logged
    Given "Data" category in context panel is expanded
    When user takes a snapshot of timelines viewer
    And user selects "AESEV" in "Split By" property in context panel
    Then "Split By" property of timelines viewer should be "AESEV"
    And timelines viewer should have repainted
    And the legend of timelines viewer should list 15 items
    And no errors should have been logged
    When user clicks on "CARDIAC DISORDERS" item in the legend of timelines viewer
    Then "CARDIAC DISORDERS" legend item in legend of timelines viewer should be selected
    And the "intervals" reading of timelines viewer should be 11
    And 143 rows should pass the filter
    And no errors should have been logged
    When user takes a snapshot of timelines viewer
    And user selects "USUBJID" in "Split By" property in context panel
    Then "Split By" property of timelines viewer should be "USUBJID"
    And timelines viewer should have repainted
    And no errors should have been logged

  Scenario: Reset View brings back every interval after a zoom
    When user remembers the "intervals" reading of timelines viewer
    And user scrolls the mouse wheel up 5 times over the "view" area of timelines viewer
    Then the "intervals" reading of timelines viewer should be lower than remembered
    When user picks "Reset View" from the context menu of timelines viewer
    Then the "intervals" reading of timelines viewer should be as remembered
    And no errors should have been logged
