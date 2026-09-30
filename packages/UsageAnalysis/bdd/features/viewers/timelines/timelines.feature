@viewers @realizes:charts.viewer.timelines
Feature: Timelines legend, legend visibility, legend clicks and Reset View
  The Timelines viewer draws one lane per subject (Split By) with an interval per event, colored by
  Color. It reports no areas or readings of its own: its legend is the standard DOM legend, so its
  items are counted and clicked by name, and what it draws is claimed by the ink of its canvas.
  ae.csv (143 adverse events, 71 subjects; AESOC has 15 values) is opened from Browse and the viewer
  is added from the ribbon's Add viewer gallery.
  "Not blank" after a legend click is claimed by the hues on its canvases: the intervals are drawn in
  the Color column's colors, and a plot area with no interval left (lane labels, axes and the zoom
  sliders only) shows 11 hue buckets against 17-18 with one category drawn, so the claim asks for 14.
  Translated from the TestTrack case Charts/timelines. Kept without (see the request document): that
  the Color list offers only string columns — no step reads the choices of a column property; that
  deselecting the last legend item draws all events again (claimed as more ink than the one-item
  state, not as the full picture); that Reset View repaints — the viewer reports no hit areas, so no
  step can zoom it first, and Reset View on an untouched view is claimed only as still painted.

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
    When user takes a snapshot of timelines viewer
    And user clicks on "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" item in the legend of timelines viewer
    Then timelines viewer should have less ink than before
    And the canvases of timelines viewer should be painted in at least 14 colors
    And 143 rows should pass the filter
    And no errors should have been logged
    When user takes a snapshot of timelines viewer
    And user clicks on "NERVOUS SYSTEM DISORDERS" item in the legend of timelines viewer holding Control
    Then timelines viewer should have more ink than before
    And no errors should have been logged
    When user takes a snapshot of timelines viewer
    And user clicks on "NERVOUS SYSTEM DISORDERS" item in the legend of timelines viewer holding Control
    Then timelines viewer should have less ink than before
    When user takes a snapshot of timelines viewer
    And user clicks on "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" item in the legend of timelines viewer
    Then timelines viewer should have more ink than before
    And "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" legend item in legend of timelines viewer should not be selected
    And no errors should have been logged
    Given "Data" category in context panel is expanded
    When user takes a snapshot of timelines viewer
    And user selects "AESEV" in "Split By" property in context panel
    Then "Split By" property of timelines viewer should be "AESEV"
    And timelines viewer should have repainted
    And the legend of timelines viewer should list 15 items
    And no errors should have been logged
    When user takes a snapshot of timelines viewer
    And user clicks on "CARDIAC DISORDERS" item in the legend of timelines viewer
    Then "CARDIAC DISORDERS" legend item in legend of timelines viewer should be selected
    And timelines viewer should have less ink than before
    And the canvases of timelines viewer should be painted in at least 14 colors
    And 143 rows should pass the filter
    And no errors should have been logged
    When user takes a snapshot of timelines viewer
    And user selects "USUBJID" in "Split By" property in context panel
    Then "Split By" property of timelines viewer should be "USUBJID"
    And timelines viewer should have repainted
    And no errors should have been logged

  Scenario: Reset View leaves the viewer painted
    When user picks "Reset View" from the context menu of timelines viewer
    Then timelines viewer should be painted
    And no errors should have been logged
