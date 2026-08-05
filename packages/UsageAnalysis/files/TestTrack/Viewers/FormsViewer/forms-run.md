# Forms viewer (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Forms from the Viewers toolbox | PASS | A field per column; one form on screen |
| 2 | The form shows the current row and follows it | PASS | Values match the table for the current row, and again after moving to row 77 |
| 3 | Show Selected Rows adds a form per selected row | PASS | Three selected rows → at least three forms; the values shown belong to one of them |
| 4 | Show Current Row hides the current-row form | PASS | 1 → 0 → 1 forms |
| 5 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 11s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* The Forms viewer's root carries **no `name` attribute**, so it is addressed by
  its class `.d4-multi-form`. The shared title-bar and gear helpers key off
  `[name="viewer-…"]`, so this spec resolves both through the root's `.panel-base`
  ancestor instead.
* Forms are counted by how many times one field name repeats on screen — the
  viewer renders one `[name="input-SEX"]` per row form.
* **Show Selected Rows ships enabled.** A blind toggle turns it off and the step
  then measures the opposite of what it claims; the shared
  `setPropertyGridCheckbox` helper was added for exactly this and is used here.

## Not automated

* Field ordering and the Fields / Sort By editors, colour coding, renderer size
  and number format.
