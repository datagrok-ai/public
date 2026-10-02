@viewers @realizes:viewers.grid
Feature: Grid summary columns
  The renderer items of the grid's Add > Summary Columns menu are this package's: every cell
  renderer it registers with `gridChart` and `virtual` set (`src/package.ts`) is offered there under
  its friendly name, and picking one appends a virtual column of that name. The grid reports the
  column in `column order` and the renderer it was given as `cell type of <col>` (the function's
  cell type, set before any drawing), so each is paired with a claim that a cell of the new column
  is painted inside its gridlines, on row 3: the menu's right click on row 2 made that row current
  and tints it, and an empty cell of a row neither current nor hovered reads blank (checked
  against an unflagged Tags cell on 2026-09-23). Also this package's: the handler that follows a
  source column into its summary columns when it is renamed or removed. Every scenario starts on a fresh
  demog-1000, whose own columns are USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG,
  CONTROL, STARTED and SEVERITY.

  Sources: `grid-appearance-summary-persist.md` scenario 6 (the one-click types) and scenario 8 (the
  summary columns across a layout and a project), `grid.md` "Removing a Summary Column by the
  Top-Panel Icon" (GROK-18256) and "Summary Column Surviving a Removed Source Column" (GROK-19942).
  The rename scenario is beyond the TestTrack specs. The menu offers eight renderer items — the
  seven of the TestTrack spec and VlaaiVis (operator decision D3: all of them); the `bar` renderer
  the package also registers is not virtual and the menu does not list it (measured on dev; not
  claimed, since it would be listed under a friendly name, not as "bar"). The three form items
  above them belong to the core and are in `UsageAnalysis/bdd/features/viewers/grid/grid-forms-column.feature`.
  The layout is saved to the server and loaded over a fresh view of the table with a scatter plot
  added; the project is saved and opened with the library's API steps (operator decision D1). A
  summary cell's own content is claimed through the Context Panel, which lists the value of every
  column the summary column draws, and the source column is renamed in its Column Properties dialog.
  The Background waits for the package autostarts: the package subscribes its rename and remove
  handler to every grid added after its own autostart, which lands seconds after the shell; a
  table opened before it got no handler until the review of 2026-09-22 (the autostart now attaches
  to the open table views too), and the wait keeps the scenarios independent of that timing. The
  GROK-19942 scenario removes its Tags column before it ends: left on the view whose CONTROL was
  removed, it made the next table opened on the page log "NullError: method not found: 'Q' on null"
  (4 runs of 4 on dev, never without that scenario before; reported to the operator). The package's
  rename and remove handler used to skip a Tags column: it looked the settings up under the grid
  column's cell type, which the grid reports in lower case, where the Tags renderer keeps them under
  `Tags` — fixed in review on 2026-09-22 (`src/package.ts` resolves the renderer's own key); the
  settings are not read here, since a Tags cell puts nothing in the Context Panel, so the
  GROK-19942 claim still rests on the renderer skipping a missing column. Which rows a Tags
  column marks is the last scenario: it carried `@known-failure` for GROK-20888 (a Tags column
  painted the grid over in one colour) until the fix of 2026-09-15 reached the stand; the marked
  cell is claimed by its hue — the chip is one colour on white — and the flood by the unmarked
  cells of row 1 staying free of any hue (row 1 is neither flagged nor current: the menu's right
  click made row 2 current). The grid
  reports `cell type of <col>` for the columns it draws, so the round-trip scenario reads the first
  three summary columns, walks the current cell to the last one so the rest scroll into view, and
  reads them; the current column is not claimed there, since a virtual column is no table column.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"

  Scenario Outline: Add > Summary Columns > <item> appends a <type> column
    When user picks "Add > Summary Columns > <item>" from the context menu of the "cell 2 of USUBJID" area of grid
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, <item>"
    And the "cell type of <item>" reading of grid should be "<type>"
    And the "cell 3 of <item>" area of grid should be painted
    And the table should have 11 columns
    And no errors should have been logged

    Examples:
      | item                | type               |
      | Sparklines          | sparkline          |
      | Bar Chart           | barchart           |
      | Pie Chart           | piechart           |
      | VlaaiVis            | vlaaivis           |
      | Radar               | radar              |
      | Smart Form          | smartform          |
      | Tags                | tags               |
      | Confidence Interval | confidenceinterval |

  Scenario: The menu lists the renderer items
    When user right-clicks on the "cell 2 of USUBJID" area of grid
    And user hovers over "Add" menu item in context menu
    And user hovers over "Summary Columns" menu item in context menu
    Then "Sparklines" menu item in context menu should be visible
    And "Bar Chart" menu item in context menu should be visible
    And "VlaaiVis" menu item in context menu should be visible
    And "Confidence Interval" menu item in context menu should be visible
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The remove icon on the top panel removes a selected summary column (GROK-18256)
    When user picks "Add > Summary Columns > Sparklines" from the context menu of the "cell 2 of USUBJID" area of grid
    Then grid should have a "header Sparklines" area
    And remove selected columns icon should be disabled
    When user clicks on the "header Sparklines" area of grid holding Control
    Then remove selected columns icon should be enabled
    When user clicks on remove selected columns icon
    Then grid should not have a "header Sparklines" area
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And the table should have 11 columns
    And no errors should have been logged

  Scenario: The grid keeps drawing a Tags column whose source column was removed (GROK-19942)
    When user picks "Add > Summary Columns > Tags" from the context menu of the "cell 2 of USUBJID" area of grid
    Then the "cell type of Tags" reading of grid should be "tags"
    And the "cell 3 of Tags" area of grid should be painted in at least 1 colors
    When user picks "Remove" from the context menu of the "header CONTROL" area of grid
    Then the table should not have a column "CONTROL"
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, STARTED, SEVERITY, Tags"
    When user clicks on the "cell 3 of AGE" area of grid
    And user presses Control+End
    Then grid should have repainted
    And grid should have a "cell 1000 of Tags" area
    And the "cell type of Tags" reading of grid should be "tags"
    When user presses Control+Home
    Then grid should have a "cell 1 of Tags" area
    And no errors should have been logged
    When user picks "Remove" from the context menu of the "header Tags" area of grid
    Then grid should not have a "header Tags" area

  Scenario: A renamed source column stays in its summary column
    When user picks "Add > Summary Columns > Sparklines" from the context menu of the "cell 2 of USUBJID" area of grid
    And user clicks on the "cell 1 of Sparklines" area of grid
    Then context panel should contain text "HEIGHT: 174.705"
    When user picks "Column Properties..." from the context menu of the "header HEIGHT" area of grid
    And user types "STATURE" into "New name:" input in HEIGHT dialog
    And user clicks on OK button in HEIGHT dialog
    Then the table should have a column "STATURE"
    And the table should not have a column "HEIGHT"
    When user clicks on the "cell 2 of Sparklines" area of grid
    Then context panel should contain text "STATURE: 150.288"
    And context panel should contain text "AGE: 30"
    And context panel should contain text "WEIGHT: 64"
    And no errors should have been logged

  Scenario: The eight summary columns and their renderers come back from a layout and a project
    When user picks "Add > Summary Columns > Sparklines" from the context menu of the "cell 2 of USUBJID" area of grid
    And user picks "Add > Summary Columns > Bar Chart" from the context menu of the "cell 2 of Sparklines" area of grid
    And user picks "Add > Summary Columns > Pie Chart" from the context menu of the "cell 2 of Bar Chart" area of grid
    And user picks "Add > Summary Columns > VlaaiVis" from the context menu of the "cell 2 of Pie Chart" area of grid
    And user picks "Add > Summary Columns > Radar" from the context menu of the "cell 2 of VlaaiVis" area of grid
    And user picks "Add > Summary Columns > Smart Form" from the context menu of the "cell 2 of Radar" area of grid
    And user picks "Add > Summary Columns > Tags" from the context menu of the "cell 2 of Smart Form" area of grid
    And user picks "Add > Summary Columns > Confidence Interval" from the context menu of the "cell 2 of Tags" area of grid
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines, Bar Chart, Pie Chart, VlaaiVis, Radar, Smart Form, Tags, Confidence Interval"
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    Then grid should not have a "header Sparklines" area
    When user loads the saved layout
    Then scatter plot viewer should be absent
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines, Bar Chart, Pie Chart, VlaaiVis, Radar, Smart Form, Tags, Confidence Interval"
    And the "cell type of Sparklines" reading of grid should be "sparkline"
    And the "cell type of Bar Chart" reading of grid should be "barchart"
    And the "cell type of Pie Chart" reading of grid should be "piechart"
    When user clicks on the "cell 2 of Pie Chart" area of grid
    And user presses ArrowRight
    And user presses ArrowRight
    And user presses ArrowRight
    And user presses ArrowRight
    And user presses ArrowRight
    Then grid should have a "header Confidence Interval" area
    And the "cell type of VlaaiVis" reading of grid should be "vlaaivis"
    And the "cell type of Radar" reading of grid should be "radar"
    And the "cell type of Smart Form" reading of grid should be "smartform"
    And the "cell type of Tags" reading of grid should be "tags"
    And the "cell type of Confidence Interval" reading of grid should be "confidenceinterval"
    When user saves the current view as project "zz-grid-summary-columns"
    And user closes all views
    And user opens the "zz-grid-summary-columns" project
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines, Bar Chart, Pie Chart, VlaaiVis, Radar, Smart Form, Tags, Confidence Interval"
    And the "cell type of Sparklines" reading of grid should be "sparkline"
    And the "cell type of Bar Chart" reading of grid should be "barchart"
    And the "cell type of Pie Chart" reading of grid should be "piechart"
    When user clicks on the "cell 2 of Pie Chart" area of grid
    And user presses ArrowRight
    And user presses ArrowRight
    And user presses ArrowRight
    And user presses ArrowRight
    And user presses ArrowRight
    Then grid should have a "header Confidence Interval" area
    And the "cell type of VlaaiVis" reading of grid should be "vlaaivis"
    And the "cell type of Radar" reading of grid should be "radar"
    And the "cell type of Smart Form" reading of grid should be "smartform"
    And the "cell type of Tags" reading of grid should be "tags"
    And the "cell type of Confidence Interval" reading of grid should be "confidenceinterval"
    And no errors should have been logged

  Scenario: A Tags column marks the rows that carry the flag and leaves the rest of the grid alone (GROK-20888)
    When user picks "Add > Summary Columns > Sparklines" from the context menu of the "cell 2 of USUBJID" area of grid
    And user picks "Add > Summary Columns > Tags" from the context menu of the "cell 2 of Sparklines" area of grid
    Then the "cell type of Tags" reading of grid should be "tags"
    And the "cell 1 of Tags" area of grid should be painted in no color
    And the "cell 3 of Tags" area of grid should be painted in at least 1 colors
    And the "cell 1 of Tags" and "cell 3 of Tags" areas of grid should be painted in different colors
    And the "cell 1 of SEVERITY" area of grid should be painted in no color
    And the "text of cell 1 of SEVERITY" reading of grid should be "High"
    And the "text of cell 1 of CONTROL" reading of grid should be "false"
    And the "text of cell 3 of CONTROL" reading of grid should be "true"
    And no errors should have been logged
