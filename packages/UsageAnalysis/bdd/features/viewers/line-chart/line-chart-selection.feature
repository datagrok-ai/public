@journey @viewers @realizes:viewers.line-chart
Feature: Line chart selection and the row markers
  What a Shift-drag over the plot selects, how the selection is painted, and the four flags that
  decide which rows the chart marks.
  The old spec counted "selection hue" pixels — pixels whose hue was near the selection colour —
  and compared the counts before and after; the chart reports `rows selected` and the library reads
  the highlight with a margin, so the claim is the rows and the paint together. To dismiss the
  context menu it opened, the old spec scanned the top row of the window at 40-pixel steps for a
  point outside every popup and over nothing clickable, clicked there, and asserted the selection
  had not changed — that scan is the platform's own "close the menu", and the claim it stood for
  is a scenario of its own here.
  **The old "Lasso selection" scenario is not translated, because the line chart has no lasso.**
  `_initAreaSelection` (line_chart_core.dart:744) makes a Shift-drag an X-BAND selection: it takes
  `xToWorld(bounds.left)` and `xToWorld(bounds.right)` and selects every row whose X falls between
  them, whatever the drag's vertical extent or shape. `lassoTool` is the annotation-region drawing
  mode ("lasso region drawing mode instead of polygon", line_chart_look.dart:430), not a selection
  tool — the old scenario turned it on, drew a three-segment zigzag, and read the band between the
  zigzag's first and last x as proof that "the lasso selects the points inside a polygon". Measured
  here: the same Shift-drag selects the same 82 rows with `lassoTool` on and off, and the library's
  own lasso gesture — which returns to its starting point — selects 0, because its band is empty.
  Fixture: spgi-100 on `CAST Idea ID` × `Chemical Space X`, one marker per row.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName      | CAST Idea ID     |
      | yColumnNames     | Chemical Space X |
      | showSelectedRows | true             |
    Then 100 rows should pass the filter
    And the "markers drawn" reading of line chart viewer should be 100
    And the "rows selected" reading of line chart viewer should be 0
    And line chart viewer should show no selection highlight
    And line chart viewer should report no error

  Scenario: A Shift-drag selects every row whose X falls in the band it spans
    Given user clears the row selection
    When user drags a selection box over the "plot" area of line chart viewer
    Then the "rows selected" reading of line chart viewer should be 82
    And 82 rows should be selected
    And line chart viewer should show more selection highlight than before
    And line chart viewer should have repainted
    When user drags a deselection box over the "plot" area of line chart viewer
    Then the "rows selected" reading of line chart viewer should be 0
    And line chart viewer should show no selection highlight
    And no errors should have been logged

  Scenario: The Lasso Tool property does not change what a Shift-drag selects
    Given user clears the row selection
    When user sets "lassoTool" property of line chart viewer to "true"
    And user drags a selection box over the "plot" area of line chart viewer
    Then the "rows selected" reading of line chart viewer should be 82
    When user clears the row selection
    And user sets "lassoTool" property of line chart viewer to "false"
    And user drags a selection box over the "plot" area of line chart viewer
    Then the "rows selected" reading of line chart viewer should be 82
    When user clears the row selection
    Then no errors should have been logged

  Scenario: Lasso Tool is offered under Tools, and picking it flips the property
    Then "lassoTool" property of line chart viewer should be "false"
    When user picks "Tools > Lasso Tool" from the context menu of line chart viewer
    Then "lassoTool" property of line chart viewer should be "true"
    When user picks "Tools > Lasso Tool" from the context menu of line chart viewer
    Then "lassoTool" property of line chart viewer should be "false"
    And no errors should have been logged

  Scenario: A drag that ends where it started selects an empty band
    Given user clears the row selection
    When user drags a lasso over the "plot" area of line chart viewer
    Then the "rows selected" reading of line chart viewer should be 0
    And line chart viewer should show no selection highlight
    And no errors should have been logged

  Scenario: Show Selected Rows decides whether the selection is painted at all
    When user selects rows where "Stereo Category" is "R_ONE"
    Then the "rows selected" reading of line chart viewer should be 36
    And line chart viewer should show a selection highlight
    When user sets "showSelectedRows" property of line chart viewer to "false"
    Then line chart viewer should show less selection highlight than before
    And the "rows selected" reading of line chart viewer should be 36
    When user sets "showSelectedRows" property of line chart viewer to "true"
    Then line chart viewer should show more selection highlight than before
    When user clears the row selection
    Then the "rows selected" reading of line chart viewer should be 0
    And no errors should have been logged

  Scenario: Opening and dismissing the context menu leaves the selection and the current row alone
    When user selects rows where "Stereo Category" is "S_ABS"
    And user makes row 7 current
    Then the "rows selected" reading of line chart viewer should be 2
    When user opens the context menu of line chart viewer
    And user closes the context menu
    Then the "rows selected" reading of line chart viewer should be 2
    And row 7 should be current
    And 2 rows should be selected
    When user clears the row selection
    Then no errors should have been logged

  Scenario: The current-row line is drawn only when it is asked for
    When user makes row 40 current
    And user sets "showCurrentRowLine" property of line chart viewer to "true"
    Then line chart viewer should have more ink than before
    When user makes row 80 current
    Then line chart viewer should have repainted
    When user sets "showCurrentRowLine" property of line chart viewer to "false"
    Then line chart viewer should have less ink than before
    And no errors should have been logged

  Scenario: Selecting through the data and reading it off the chart agree
    When user selects rows where "Stereo Category" is one of "S_ABS, S_PART"
    Then 12 rows should be selected
    And the "rows selected" reading of line chart viewer should be 12
    And line chart viewer should show a selection highlight
    When user selects no rows
    Then the "rows selected" reading of line chart viewer should be 0
    And line chart viewer should show no selection highlight
    And no errors should have been logged

  @known-failure
  Scenario: After Tools > Lasso Tool the next Shift-drag selects nothing
    # OPEN BUG, no ticket yet. Setting `lassoTool` through the property leaves the Shift-drag
    # selecting its usual 82 rows (the scenario above this one); picking the SAME property from
    # the chart's own Tools menu leaves the next Shift-drag selecting 0. It is not the drawing
    # mode — right after the pick the chart reports `region drawing mode` false and
    # `viewer regions` 0, so `_initAreaSelection`'s `isInDrawingMode` guard is not the one that
    # is closing (line_chart_core.dart:747). Reproduced on every run; the API path and the UI
    # path for one property must not differ.
    Given user clears the row selection
    When user picks "Tools > Lasso Tool" from the context menu of line chart viewer
    Then "lassoTool" property of line chart viewer should be "true"
    And the "region drawing mode" reading of line chart viewer should be "false"
    When user drags a selection box over the "plot" area of line chart viewer
    Then the "rows selected" reading of line chart viewer should be 82
    When user clears the row selection
    And user sets "lassoTool" property of line chart viewer to "false"
    Then no errors should have been logged
