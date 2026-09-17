@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot selection, current row and mouse-over
  The mouse on the chart: a Shift-drag over the band between two axes selects every polyline it
  crosses and is additive by design (no Control needed), a click on chart space no polyline passes
  through clears the selection, a click on a polyline makes that row current, and a hover makes it
  the mouse-over row. Show All Lines off then leaves only the current, hovered and selected
  polylines, which is what `lines drawn` counts. Every gesture aims at a region the plot reports —
  `band <colA> - <colB>`, `line of row <n>`, `empty space` — so nothing here clicks a fraction of
  the canvas.
  demog-1000 with AGE, HEIGHT and WEIGHT: 1000 rows, 872 of which have a HEIGHT, 896 Caucasian and
  15 Asian. The band scenarios first filter the blank HEIGHTs out, because the plot's rectangle
  selection reads a blank as the raw storage sentinel rather than skipping it (see the defect note
  below), so on the full table a Shift-drag also catches rows whose line is not drawn between those
  two axes. With the blanks filtered, a drag over the inner 80% of the AGE-HEIGHT band selects
  exactly the 872 rows on screen. Two of the 15 Asian rows (899 and 902) have no HEIGHT, so they are
  outside that filter and survive the drag only because the drag adds.
  Opening the table makes row 1 current, so `current row` reads 1 and `lines drawn` never falls
  below that one line. Every scenario clears what it selected and puts the current row back.
  Defect found while writing this: `_belongsToRect` and `_hitTest` call `col.toDouble(row)` without
  the `isFinite` guard `_renderLines` uses, so a row with a blank value is hit-tested at ~0 instead
  of being skipped — on the full table 61 of the 128 blank-HEIGHT rows are selected by a drag that
  never crossed their line.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then pc plot viewer should show 1000 rows
    And the "rows selected" reading of pc plot viewer should be 0
    And the "current row" reading of pc plot viewer should be 1
    And table "demog-1000" should have missing values in "HEIGHT" column
    And pc plot viewer should have a "band \"AGE\" - \"HEIGHT\"" area
    And pc plot viewer should have a "line of row 1000" area

  Scenario: A Shift-drag over a band selects the polylines it crosses
    When user filters rows where "HEIGHT" is not null
    Then 872 rows should pass the filter
    And pc plot viewer should show 872 rows
    When user drags a selection box over the "band \"AGE\" - \"HEIGHT\"" area of pc plot viewer
    Then 872 rows should be selected
    And the "rows selected" reading of pc plot viewer should be 872
    And the line of row 1000 of pc plot viewer should be selected
    And every selected row should pass the filter
    And pc plot viewer should show more selection highlight than before
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: A second Shift-drag adds to the selection instead of replacing it
    Given 872 rows should pass the filter
    When user selects rows where "RACE" is "Asian"
    Then 15 rows should be selected
    When user drags a selection box over the "band \"AGE\" - \"HEIGHT\"" area of pc plot viewer
    Then 874 rows should be selected
    And all rows where "RACE" is "Asian" should be selected
    And the line of row 899 of pc plot viewer should be selected
    And the "rows selected" reading of pc plot viewer should be higher than before
    When user clears the row selection
    And user resets the filter
    Then all rows should pass the filter
    And no rows should be selected
    And pc plot viewer should show 1000 rows
    And no errors should have been logged

  Scenario: A click on chart space no polyline passes through clears the selection
    When user selects rows where "RACE" is "Caucasian"
    Then 896 rows should be selected
    And pc plot viewer should show more selection highlight than before
    When user clicks on the "empty space" area of pc plot viewer
    Then no rows should be selected
    And the "rows selected" reading of pc plot viewer should be 0
    And pc plot viewer should show less selection highlight than before
    And no errors should have been logged

  Scenario: A click on a polyline makes that row current
    Given user listens for "d4-pc-plot-on-line-clicked" event on pc plot viewer
    And no rows should be selected
    Then the "current row" reading of pc plot viewer should be 1
    When user clicks on the "line of row 1000" area of pc plot viewer
    Then "d4-pc-plot-on-line-clicked" event should have fired on pc plot viewer
    And the "current row" reading of pc plot viewer should be 1000
    And row 1000 should be current
    And no rows should be selected
    When user clicks on the "line of row 999" area of pc plot viewer
    Then the "current row" reading of pc plot viewer should be 999
    When user makes row 1 current
    Then the "current row" reading of pc plot viewer should be 1
    And no errors should have been logged

  Scenario: Hovering a polyline makes it the mouse-over row
    Given user listens for "d4-pc-plot-on-line-hovered" event on pc plot viewer
    When user hovers over the "empty space" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 0
    When user hovers over the "line of row 1000" area of pc plot viewer
    Then "d4-pc-plot-on-line-hovered" event should have fired on pc plot viewer
    And the "hovered row" reading of pc plot viewer should be 1000
    When user hovers over the "empty space" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 0
    And no errors should have been logged

  Scenario: Show All Lines off leaves the current, selected and hovered polylines
    When user clears the row selection
    And user hovers over the "empty space" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 0
    And the "lines drawn" reading of pc plot viewer should be 1000
    When user sets "Show All Lines" property of pc plot viewer to "false"
    Then the "lines drawn" reading of pc plot viewer should be 1
    And pc plot viewer should have less ink than before
    When user selects rows where "RACE" is "Asian"
    Then the "lines drawn" reading of pc plot viewer should be 16
    And pc plot viewer should have more ink than before
    When user hovers over the "line of row 1000" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 1000
    And the "lines drawn" reading of pc plot viewer should be 17
    When user clears the row selection
    Then the "lines drawn" reading of pc plot viewer should be 2
    When user hovers over the "empty space" area of pc plot viewer
    And user sets "Show All Lines" property of pc plot viewer to "true"
    Then the "lines drawn" reading of pc plot viewer should be 1000
    And pc plot viewer should have more ink than before
    And no errors should have been logged

  Scenario: The Selection menu drives the lines and the overlay
    When user clears the row selection
    And user hovers over the "empty space" area of pc plot viewer
    Then the "lines drawn" reading of pc plot viewer should be 1000
    When user picks "Selection > Show All Lines" from the context menu of pc plot viewer
    Then "Show All Lines" property of pc plot viewer should be "false"
    And pc plot viewer should have less ink than before
    When user picks "Selection > Show Current Line" from the context menu of pc plot viewer
    Then "Show Current Line" property of pc plot viewer should be "false"
    And pc plot viewer should have repainted by at least 200 pixels
    When user picks "Selection > Show Current Line" from the context menu of pc plot viewer
    Then "Show Current Line" property of pc plot viewer should be "true"
    And pc plot viewer should have repainted by at least 200 pixels
    When user picks "Selection > Show Mouse Over Line" from the context menu of pc plot viewer
    Then "Show Mouse Over Line" property of pc plot viewer should be "false"
    When user picks "Selection > Show Mouse Over Line" from the context menu of pc plot viewer
    Then "Show Mouse Over Line" property of pc plot viewer should be "true"
    When user picks "Selection > Show All Lines" from the context menu of pc plot viewer
    Then "Show All Lines" property of pc plot viewer should be "true"
    And the "lines drawn" reading of pc plot viewer should be 1000
    And pc plot viewer should have more ink than before
    And no errors should have been logged
