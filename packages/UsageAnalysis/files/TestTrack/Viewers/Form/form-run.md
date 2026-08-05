# Form viewer (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Form from the Viewers toolbox | PASS | A field is rendered for every column of the table |
| 2 | The fields show the current row | PASS | USUBJID / SEX / RACE / DIS_POP match the table verbatim |
| 3 | The next and previous arrows walk the rows | PASS | Current row moves ±1 and the fields follow |
| 4 | Picking a row in the grid updates the form | PASS | Row 42, fields match |
| 5 | The row selector selects the row on show | PASS | Exactly one row selected, then cleared |
| 6 | Show Navigation hides the toolbar | PASS | Toolbar hidden and restored |
| 7 | Show Next Row Arrow hides just that arrow | PASS | Next hidden while previous stays visible |
| 8 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 11s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* This viewer is plain DOM: fields are `[name="div-<COLUMN>"]` labels with
  `[name="input-<COLUMN>"]` values (underscores in a column name become hyphens).
  That makes the strongest assertion available on any of these viewers possible —
  the text on screen is compared verbatim with the table's own values.
* Only string columns are compared literally; numbers and dates are shown
  formatted (`160.484`, `9/11/1990`), so those fields are covered by the
  "it changed with the row" checks instead.

## Not automated

* The field editor, the built-in form designer (drag-and-drop layout), and
  save/open of a form file — `form-ui.md`.
