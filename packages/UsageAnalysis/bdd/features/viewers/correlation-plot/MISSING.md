# Correlation plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `correlation_plot/CLAUDE.md` "Automation
surface" and the TestTrack `correlation-plot.md` / `correlation-plot-ui.md`; trimmed in the BDD PR round of 8 Nov 2026.

## Resolved

- Selection highlighting in the diagonal histograms (core#1087): a scenario in `correlation-plot-colors-and-chrome.feature`.
- A column removed from narrowed axes leaves only that column (core#1172): Dart `Correlation plot | Removing a column
  takes only that column off the axes`.

## Still open

- **Grid > Grid Color Coding > Apply to text** (GROK-19052): the coloured text is drawn on the inner grid's canvas;
  needs a `text color of cell <X> x <Y>` reading (`GridCellStyle.textColor` at render time, in `correlation_status.dart`).

## Won't do

- No row passing, Row Source Selected with one row, a table with no numeric column, 120 numeric columns opening without
  coefficients, a column renamed on the axes: not UI checks (Dart `Data Manipulation`; renames: ApiTests
  `AI: Viewers: Shared rules`).
