@journey @viewers @realizes:viewers.line-chart
Feature: Line chart statistical process control, zoom and Reset View
  The SPC overlay — the control limits, the two sigma bands, the process average and the Western
  Electric rule violations — and what the wheel and Reset View do to the X window.
  The old SPC scenario's entire claim that the chart survived being switched on was
  `await page.evaluate(() => true)`, labelled "page stays responsive, no freeze (GROK-20126)". A
  round trip to the page proves the page exists; it says nothing about whether SPC drew. The chart
  now reports `upper control limit`, `lower control limit`, `spc average` and `violations`, and the
  bands are hit areas that are absent unless the SPC pass drew them — so hiding a band is the band
  disappearing, and a hand-set limit is the violation count moving.
  The zoom half was half honest already (it paced five WheelEvents and subscribed to the events in
  the page) but never showed that the zoom changed anything or that the reset undid it — it called
  `resetView()` programmatically and then asserted the two X bound *properties* were untouched.
  Here the wheel is a real wheel over the plot, Reset View is the menu item, and the claim is the
  X window.
  Fixture: spgi-100 on `CAST Idea ID` × `Chemical Space X`. SPC is gated on a single un-split
  series, so the split scenario is the gate itself.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then 100 rows should pass the filter
    And "showStatisticalProcessControl" property of line chart viewer should be "false"
    And line chart viewer should not have a "control limits" area
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    And line chart viewer should report no error

  Scenario: Switching SPC on draws the limits, the sigma bands and the average
    When user sets "showStatisticalProcessControl" property of line chart viewer to "true"
    Then line chart viewer should have a "control limits" area
    And line chart viewer should have a "sigma 1" area
    And line chart viewer should have a "sigma 2" area
    And line chart viewer should have an "average" area
    And the "spc average" reading of line chart viewer should be between 3.4 and 3.6
    And the "upper control limit" reading of line chart viewer should be between 21 and 22
    And the "lower control limit" reading of line chart viewer should be between -15 and -14
    And the "violations" reading of line chart viewer should be 0
    And the "control limits" area of line chart viewer should be painted
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be higher than before
    And line chart viewer should have repainted
    When user sets "showStatisticalProcessControl" property of line chart viewer to "false"
    Then line chart viewer should not have a "control limits" area
    And line chart viewer should not have an "average" area
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be lower than before
    And line chart viewer should have repainted
    And no errors should have been logged

  Scenario: SPC is gated on a single un-split series
    When user sets "showStatisticalProcessControl" property of line chart viewer to "true"
    Then line chart viewer should have a "control limits" area
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then line chart viewer should not have a "control limits" area
    And line chart viewer should not have a "sigma 1" area
    And line chart viewer should not have an "average" area
    And the "lines" reading of line chart viewer should be 5
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then line chart viewer should have a "control limits" area
    When user sets "multiAxis" property of line chart viewer to "true"
    Then line chart viewer should not have a "control limits" area
    When user sets properties of line chart viewer:
      | multiAxis                     | false |
      | showStatisticalProcessControl | false |
    Then no errors should have been logged

  Scenario: The bands can be hidden one at a time and the limits stay reported
    When user sets "showStatisticalProcessControl" property of line chart viewer to "true"
    Then line chart viewer should have a "sigma 1" area
    When user sets "showSigma1" property of line chart viewer to "false"
    Then line chart viewer should not have a "sigma 1" area
    And line chart viewer should have a "sigma 2" area
    And the "plot" area of line chart viewer should have less ink than before
    When user sets "showSigma2" property of line chart viewer to "false"
    Then line chart viewer should not have a "sigma 2" area
    When user sets "showAverage" property of line chart viewer to "false"
    Then line chart viewer should not have an "average" area
    And the "spc average" reading of line chart viewer should be between 3.4 and 3.6
    When user sets "showControlLimits" property of line chart viewer to "false"
    Then line chart viewer should not have a "control limits" area
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be lower than before
    When user sets properties of line chart viewer:
      | showSigma1        | true |
      | showSigma2        | true |
      | showAverage       | true |
      | showControlLimits | true |
    Then line chart viewer should have a "control limits" area
    And the 'y axis max of "Chemical Space X"' reading of line chart viewer should be higher than before
    And the "control limits" area of line chart viewer should be painted
    When user sets "showStatisticalProcessControl" property of line chart viewer to "false"
    Then no errors should have been logged

  Scenario: The Western Electric rules flag points, and hand-set limits flag far more of them
    When user sets properties of line chart viewer:
      | showStatisticalProcessControl | true |
      | showBias                      | true |
      | showConsistentTrend           | true |
      | showOscillation               | true |
      | showMediumShift               | true |
      | showSustainedShift            | true |
      | showSuppressedVariation       | true |
    Then the "violations" reading of line chart viewer should be 33
    And line chart viewer should have a "violation 2" area
    And line chart viewer should have a "violation 6" area
    And line chart viewer should not have a "violation 1" area
    When user sets properties of line chart viewer:
      | lowerControlLimit | 0 |
      | upperControlLimit | 5 |
    Then the "upper control limit" reading of line chart viewer should be 5
    And the "lower control limit" reading of line chart viewer should be 0
    And the "violations" reading of line chart viewer should be 89
    And line chart viewer should have a "violation 1" area
    And the "control limits" area of line chart viewer should be shorter than before
    When user sets properties of line chart viewer:
      | lowerControlLimit |  |
      | upperControlLimit |  |
    Then the "upper control limit" reading of line chart viewer should be between 21 and 22
    And the "violations" reading of line chart viewer should be 33
    When user sets properties of line chart viewer:
      | showBias                      | false |
      | showConsistentTrend           | false |
      | showOscillation               | false |
      | showMediumShift               | false |
      | showSustainedShift            | false |
      | showSuppressedVariation       | false |
      | showStatisticalProcessControl | false |
    Then no errors should have been logged

  Scenario: The wheel zooms the X axis and Reset View puts it back
    Given user listens for "d4-linechart-zoomed" event on line chart viewer
    And user listens for "d4-linechart-reset-view" event on line chart viewer
    When user scrolls the mouse wheel up over the "plot" area of line chart viewer
    Then "d4-linechart-zoomed" event should have fired on line chart viewer
    And the "x axis span" reading of line chart viewer should be lower than before
    And line chart viewer should have repainted
    When user picks "Reset View" from the context menu of line chart viewer
    Then "d4-linechart-reset-view" event should have fired on line chart viewer
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    And no errors should have been logged

  Scenario: X Min and X Max pin the window, and Reset View returns to them rather than to the column
    When user sets properties of line chart viewer:
      | xMin | 634800 |
      | xMax | 634850 |
    Then the "x axis min" reading of line chart viewer should be 634800
    And the "x axis max" reading of line chart viewer should be 634850
    And the "x axis span" reading of line chart viewer should be 50
    And the "markers drawn" reading of line chart viewer should be lower than before
    When user scrolls the mouse wheel up over the "plot" area of line chart viewer
    Then the "x axis span" reading of line chart viewer should be lower than before
    When user picks "Reset View" from the context menu of line chart viewer
    Then "xMin" property of line chart viewer should be "634800"
    And "xMax" property of line chart viewer should be "634850"
    And the "x axis span" reading of line chart viewer should be 50
    When user sets properties of line chart viewer:
      | xMin |  |
      | xMax |  |
    Then the "x axis span" reading of line chart viewer should be between 106 and 107
    And no errors should have been logged
