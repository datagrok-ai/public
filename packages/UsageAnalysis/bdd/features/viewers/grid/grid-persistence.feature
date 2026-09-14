@journey @viewers @realizes:viewers.grid
Feature: Grid appearance and geometry across a layout and a project
  One journey that dresses the grid up and then asks two round-trips to give all of it back: four
  colour codings (Linear on AGE, Conditional on HEIGHT, Categorical on SEX, WEIGHT linked to SEX),
  the row height, the missing-value colour, the min and max stats rows, a column moved, a column
  hidden, a column widened, a column and two rows pinned and a sort. Every claim reads the grid's
  own status where it has one: the renderer-resolved `color of cell <r> of <c>`, `column order`,
  `column width of AGE`, `row height` (with the height of a drawn cell next to it), `pinned rows`,
  `sort column` and `sort direction`; the colour codings' tags and Frozen Columns are read from the
  column and the property. The three pinned rows (2, 4 and 298, unique by USUBJID) stay on top
  whatever the sort, so their cells are the ones the colour claims read after each round-trip: row
  2 is F with HEIGHT 150.288 (the blue range), row 4 is M with HEIGHT 183.83 (the red one), row 298
  is F with no HEIGHT.
  The layout is saved to the server and loaded over a fresh view of the table with a scatter plot
  added, and the fresh view is shown to differ first, so nothing the load restores was already
  there; the project is saved and opened with the library's API steps (operator decision D1). The
  summary columns of the TestTrack round-trip (its GROK-19769 guard) are PowerGrid's, in
  `packages/PowerGrid/bdd/features/grid/summary-columns.feature`.

  Whether each colouring, the row height, the missing-value colour, a header drag, a resize, a pin
  and a sort work at all is claimed by `grid-appearance.feature`, `grid-columns.feature` and
  `grid-pinning.feature`; here they set up the state, and the only claims made before the
  round-trips are the ones that make the round-trip claims able to fail (the pinned cells are
  coloured, and coloured differently). New here: WEIGHT linked to SEX through the header menu and
  the colour-coding dialog (`color-coding.feature` links through the API), hiding a column in the
  Order or Hide Columns dialog, and the two round-trips.

  What the TestTrack specs drive through the UI and this feature does not: the Row Height and the
  Missing Value Color are set as grid properties, not in the gear panel; the layout is saved with
  `dapi.layouts.save`, not from the ribbon; the project goes through the API (D1).

  Not translated, and why: whether the min and max stats rows are drawn, and what they show
  (`grid-ui.md` "Column Stats - Visual Verification") — the grid reports no reading or area for its
  special rows, so the stats-row scenario claims only that the columns stay as they were and
  nothing is logged - a floor, not a claim that could fail if the rows were missing. A project holding a table of extracted rows (`grid.md` "Extracted Rows in a
  Saved Project", GROK-19717): the library's project step saves the current view and its table,
  and a project with every open table is deferred by operator decision D2. The column hidden is RACE, not WEIGHT as in
  the TestTrack spec: WEIGHT carries the linked colouring, and a hidden column has no cells whose
  colour a round-trip could be checked by. RACE is unticked in the Order or Hide Columns dialog,
  whose list is a grid of its own: it reports its rows as `cell <r> of __name` (the column) and
  `cell <r> of x` (the tick), in the table's column order.
  The Linked source column is picked in the colour-coding dialog; the Conditional ranges are
  written through the API, as in `grid-appearance.feature`, because the ranges editor has no named
  inputs. Row 298 is scrolled into view with `user makes row 298 current` before it is right-clicked.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And "Frozen Columns" property of grid should be "1"

  Scenario: Four colour codings, the row height and the missing-value colour are set
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Conditional" from the context menu of the "header HEIGHT" area of grid
    And user colors "HEIGHT" column conditionally:
      | <160 | #0000FF |
      | >180 | #FF0000 |
    And user picks "Color Coding > Categorical" from the context menu of the "header SEX" area of grid
    Then "AGE" column should have tag ".color-coding-type" equal to "Linear"
    And "HEIGHT" column should have tag ".color-coding-type" equal to "Conditional"
    And "SEX" column should be color-coded categorically
    When user picks "Color Coding > Linked" from the context menu of the "header WEIGHT" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header WEIGHT" area of grid
    Then "Color-coding: WEIGHT" dialog should be visible
    When user selects "SEX" in "Source column" input in "Color-coding: WEIGHT" dialog
    And user clicks on CLOSE button in "Color-coding: WEIGHT" dialog
    Then the coloring of "WEIGHT" column should be linked to "SEX" column
    And the "color of cell 1 of WEIGHT" and "color of cell 1 of SEX" readings of grid should be the same
    And the "color of cell 4 of WEIGHT" and "color of cell 4 of SEX" readings of grid should be the same
    And the "color of cell 1 of WEIGHT" and "color of cell 4 of WEIGHT" readings of grid should differ
    When user sets "Row Height" property of grid to "40"
    And user sets "Missing Value Color" property of grid to "#FFAAAA"
    Then the "row height" reading of grid should be 40
    And no errors should have been logged

  Scenario: The min and max stats rows leave the columns as they were
    When user remembers the "column order" reading of grid
    And user picks "Add > Column Stats > min" from the context menu of the "cell 2 of AGE" area of grid
    And user picks "Add > Column Stats > max" from the context menu of the "cell 2 of AGE" area of grid
    Then the "column order" reading of grid should be as remembered
    And grid should show 1000 rows
    And no errors should have been logged

  Scenario: Moving, hiding, widening, pinning and sorting set the geometry
    When user drags the "header HEIGHT" area of grid to the "header DEMOG" area
    Then the "column order" reading of grid should differ from before
    When user picks "Order or Hide Columns..." from the context menu of the "cell 2 of AGE" area of grid
    Then the "text of cell 4 of __name" reading of Grid viewer in Order or Hide Columns dialog should be "RACE"
    And the "text of cell 4 of x" reading of Grid viewer in Order or Hide Columns dialog should be "true"
    When user clicks on the "cell 4 of x" area of Grid viewer in Order or Hide Columns dialog
    Then the "text of cell 4 of x" reading of Grid viewer in Order or Hide Columns dialog should be "false"
    When user clicks on CLOSE button in Order or Hide Columns dialog
    Then grid should not have a "header RACE" area
    And the "column order" reading of grid should not contain "RACE"
    When user drags the "column resizer AGE" area of grid by 60 pixels to the right
    Then the "column width of AGE" reading of grid should be higher than before
    When user picks "Pin > Pin Column" from the context menu of the "header SEX" area of grid
    Then "Frozen Columns" property of grid should be "2"
    When user picks "Pin > Pin Row" from the context menu of the "cell 2 of USUBJID" area of grid
    And user picks "Pin > Pin Row" from the context menu of the "cell 4 of USUBJID" area of grid
    And user makes row 298 current
    And user makes row 300 current
    And user picks "Pin > Pin Row" from the context menu of the "cell 298 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 3
    When user double-clicks on the "header AGE" area of grid
    And user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "ascending"
    And the "color of cell 298 of HEIGHT" reading of grid should be "#ffaaaa"
    And the "color of cell 4 of HEIGHT" reading of grid should be "#ff0000"
    And the "color of cell 2 of HEIGHT" reading of grid should be "#0000ff"
    And the "color of cell 4 of AGE" and "color of cell 298 of AGE" readings of grid should differ
    And the "color of cell 4 of SEX" and "color of cell 298 of SEX" readings of grid should differ
    And the "color of cell 298 of WEIGHT" and "color of cell 298 of SEX" readings of grid should be the same
    When user remembers the "column order" reading of grid
    And user remembers the "column width of AGE" reading of grid
    And user remembers the "row height" reading of grid
    And user remembers the "color of cell 4 of AGE" reading of grid
    And user remembers the "color of cell 298 of AGE" reading of grid
    And user remembers the "color of cell 4 of SEX" reading of grid
    And user remembers the "color of cell 298 of SEX" reading of grid
    And user remembers the "color of cell 4 of WEIGHT" reading of grid
    And user remembers the "color of cell 298 of WEIGHT" reading of grid
    Then no errors should have been logged

  Scenario: A layout saved to the server restores everything over a fresh view of the table
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    Then the "column order" reading of grid should not be as remembered
    And the "column width of AGE" reading of grid should not be as remembered
    And the "row height" reading of grid should not be as remembered
    And the "color of cell 4 of AGE" reading of grid should not be as remembered
    And "Frozen Columns" property of grid should be "1"
    And the "pinned rows" reading of grid should be 0
    And the "sort column" reading of grid should be ""
    When user loads the saved layout
    Then scatter plot viewer should be absent
    And the "column order" reading of grid should be as remembered
    And the "column width of AGE" reading of grid should be as remembered
    And the "row height" reading of grid should be as remembered
    And the "cell 4 of AGE" area of grid should be at least 36 pixels tall
    And "Frozen Columns" property of grid should be "2"
    And the "pinned rows" reading of grid should be 3
    And the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "ascending"
    And the "color of cell 4 of AGE" reading of grid should be as remembered
    And the "color of cell 298 of AGE" reading of grid should be as remembered
    And the "color of cell 4 of HEIGHT" reading of grid should be "#ff0000"
    And the "color of cell 298 of HEIGHT" reading of grid should be "#ffaaaa"
    And the "color of cell 2 of HEIGHT" reading of grid should be "#0000ff"
    And the "color of cell 4 of SEX" reading of grid should be as remembered
    And the "color of cell 298 of SEX" reading of grid should be as remembered
    And the "color of cell 4 of WEIGHT" reading of grid should be as remembered
    And the "color of cell 298 of WEIGHT" reading of grid should be as remembered
    And the coloring of "WEIGHT" column should be linked to "SEX" column
    And no errors should have been logged

  Scenario: A project round-trip restores everything
    When user saves the current view as project "zz-grid-persistence"
    And user closes all views
    And user opens the "zz-grid-persistence" project
    Then grid should show 1000 rows
    And the "column order" reading of grid should be as remembered
    And the "column width of AGE" reading of grid should be as remembered
    And the "row height" reading of grid should be as remembered
    And the "cell 4 of AGE" area of grid should be at least 36 pixels tall
    And "Frozen Columns" property of grid should be "2"
    And the "pinned rows" reading of grid should be 3
    And the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "ascending"
    And the "color of cell 4 of AGE" reading of grid should be as remembered
    And the "color of cell 298 of AGE" reading of grid should be as remembered
    And the "color of cell 4 of HEIGHT" reading of grid should be "#ff0000"
    And the "color of cell 298 of HEIGHT" reading of grid should be "#ffaaaa"
    And the "color of cell 2 of HEIGHT" reading of grid should be "#0000ff"
    And the "color of cell 4 of SEX" reading of grid should be as remembered
    And the "color of cell 298 of SEX" reading of grid should be as remembered
    And the "color of cell 4 of WEIGHT" reading of grid should be as remembered
    And the "color of cell 298 of WEIGHT" reading of grid should be as remembered
    And no errors should have been logged
