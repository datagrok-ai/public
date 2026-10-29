# Matrix plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `matrix_plot/CLAUDE.md` "Automation
surface" and `matrix-plot-ui.md` ("Interact with Plot Elements"). What could be said is in
`matrix-plot-filter-panel-and-axes.feature`.

## Points inside a cell — the product first, then the areas

`matrix-plot-ui.md` asks that hovering a point in an off-diagonal cell highlights its row in the
grid, that a click makes it current and that a drag selects the points in its region. On this build
none of the three does anything: `MatrixPlotCore` appends only the inner viewer's **canvas** to each
cell (`matrix_plot_core.dart`, `cellDiv.children.add(cellCanvas)`), so the inner viewer's mouse
handlers, bound to its own root and overlay, never receive an event. Measured on demog-1000: a plain
drag and a Shift-drag across `cell body HEIGHT x AGE` select 0 rows with Density plot and with
Scatter plot cells; a hover leaves `mouseOverRowIdx` at -1 and a click leaves the current row where
it was. The old spec's own note says the same ("live recon also found them non-functional via trusted
input"). Whether that is intended is a product question; if it is not, the scenarios need, besides the
fix, the inner scatter plot's marker areas re-exported per cell (`marker of row <n> in cell <X> x
<Y>`, from the inner viewer's `getWidgetStatus`, translated by the cell's offset in
`features/matrix_status.dart`):

```gherkin
Scenario: A point inside a cell hovers, clicks and drag-selects its row
  Given user sets "cellPlotType" property of matrix plot viewer to "Scatter plot"
  When user hovers over the "marker of row 11 in cell HEIGHT x AGE" area of matrix plot viewer
  Then the mouse-over row of the table should be 11
  When user clicks on the "marker of row 11 in cell HEIGHT x AGE" area of matrix plot viewer
  Then row 11 should be current
  When user drags a selection box over the "cell body HEIGHT x AGE" area of matrix plot viewer
  Then <measured> rows should be selected
  And no errors should have been logged
```

## Edge cases (`matrix-plot-edge-cases.feature`)

Written for `matrix-plot-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Steps requested**

1. A layout naming a column the table no longer has gives a blank plot with no message (confirmed on the stand through the API: `loadLayout` after HEIGHT was removed, or `xColumnNames` with `NO_SUCH`). There is no manual path a feature can take: the Layouts pane of the toolbox does not offer a layout once its table lacks one of its columns (measured: saved from the pane on `mp-edge-<time>`, HEIGHT removed, the table opened again without HEIGHT — the card is not listed). The scenario needs either a UI path that applies such a layout (a layout file opened on the table, or the pane listing it) or stays a Dart regTest.

**Decided**

1. A layout naming a missing column: the missing column is skipped and the others are drawn (the layout rule of `viewer-rules-edge-cases.feature`). Today nothing is drawn and nothing is said, a defect with no ticket. No outline row yet: the matrix plot has no reading of the columns it lays out after a layout, and its X Column Names may keep the missing name while drawing nothing; the cell count (`cells drawn`) is the reading to claim once the expected count is measured.

**API or Dart tests to write**

1. Dart regTest: `xColumnNames` naming a column the table lacks (`matrix_plot_core.dart` keeps the null from `dataFrame[cn]`) draws the other columns.

### Won't do

- 0-cell viewport leaves the render pending: not reproduced on the stand (`columns` 0, not pending, at 0.5, 3 and 12 s); the note in `matrix-plot-viewport.feature:30` no longer holds.
