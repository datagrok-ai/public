# Matrix plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `matrix_plot/CLAUDE.md` "Automation
surface" and `matrix-plot-ui.md` ("Interact with Plot Elements"); trimmed in the BDD PR round of 8 Nov 2026.

## Still open

- **Points inside a cell** (`matrix-plot-ui.md`): a hover, a click and a drag over a point of an off-diagonal cell do
  nothing, as `MatrixPlotCore` appends only the inner viewer's canvas to the cell (`matrix_plot_core.dart`). A product
  question first; then the inner scatter plot's marker areas re-exported per cell (`marker of row <n> in cell <X> x
  <Y>`, `features/matrix_status.dart`).
- **A layout naming a column the table no longer has draws nothing** (`matrix_plot_core.dart` keeps the null of
  `dataFrame[cn]` in the X columns) instead of the other columns — needs a ticket; claim it then in a Dart `Invalid
  settings` case.

## Won't do

- X Column Names cleared, a constant column, a column renamed, every row filtered out or deleted: the shared viewer
  rules (Dart `Data Manipulation | Matrix plot`, ApiTests `AI: Viewers: Shared rules`).
