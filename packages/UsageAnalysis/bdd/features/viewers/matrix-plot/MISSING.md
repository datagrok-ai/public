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
