# Word Cloud (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Word cloud from the Viewers toolbox | PASS | Every drawn word's count matches the column's row counts |
| 2 | Column RACE draws one word per race with its row count | PASS | 4 words, counts verified against the dataframe |
| 3 | A high-cardinality column shows an error instead of a cloud | PASS | USUBJID → "500 or fewer unique categories"; back to RACE clears it |
| 4 | Shape rearranges the words | PASS | star / diamond / circle each repaint |
| 5 | Equal min and max text size flattens the word sizes | PASS | 20/20, then restored to 12/48 |
| 6 | Rotation range changes the word angles | PASS | 0/0 repaints the cloud |
| 7 | Bold and font family restyle the words | PASS | Bold flips (ships **on**), monospace repaints |
| 8 | Hovering a word shows its row count | PASS | Tooltip matches `/\d+ rows/` |
| 9 | Clicking a word selects exactly that word's rows | PASS | Selected count equals one of the drawn word counts |
| 10 | Filtering the table leaves the cloud rendering | PASS | See the observation below |
| 11 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 38s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* The canvas is a plain 2d ECharts canvas, so repaints are proved with pixel
  histograms. What was drawn is read back from the chart's series — that is the
  same list of words and counts the user sees, and it makes the column-assignment
  steps assert real numbers instead of a property round-trip.
* Words carry no DOM nodes; hover and click targets are found as dark pixel blobs
  on the canvas.

## Observation on dev (not filed)

**The Word cloud ignores the table filter.** With `SEX = M` applied (2607 of 5850
rows), the RACE cloud keeps showing whole-table counts:

| Word | Drawn | Rows actually passing the filter |
|---|---:|---:|
| Caucasian | 5267 | 2444 |
| Other | 354 | 75 |
| Black | 157 | 53 |
| Asian | 72 | 35 |

Unchanged after 12s and after forcing a redraw with a look property. The viewer
exposes no `rowSource` and no `filter` property at all
(`columnColumnName, shape, minTextSize, maxTextSize, minRotationDegree,
maxRotationDegree, rotationStep, gridSize, drawOutOfBound, fontFamily, bold`), so
filtering appears simply not to be implemented — the same class of gap as
GROK-15543 (Sunburst: does not support collaborative filtering).

`word-cloud.md` step "Filter interaction" describes the counts updating, which
does not happen today.

## Not automated

* Grid Size, Draw Out Of Bound, viewer resize / full screen, context menu — kept
  as a manual list in `word-cloud.md`.
