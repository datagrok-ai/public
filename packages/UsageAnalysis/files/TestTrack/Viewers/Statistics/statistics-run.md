# Statistics viewer (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add, close and re-add the viewer | PASS | Content canvas populated both times |
| 2 | Statistics > sum adds a column to the viewer | PASS | Entry read as unticked, clicked, viewer redraws |
| 3 | Row Source Filtered follows the table filter | PASS | `SEX = M`, statistics redraw |
| 4 | Row Source Selected follows the selection | PASS | Every second row selected, statistics redraw |
| 5 | The Columns property lists every column of the table | PASS | `11 / 11` matches the table's column count |
| 6 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 18s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* The viewer renders its table into a canvas and its inner grid is not reachable
  from the JS API (`DG.Viewer.fromRoot` returns `Unknown`, and the viewer object
  exposes only the source dataframe). So the assertions are: the canvas redraws,
  the context-menu entry's tick state, and the property-grid values.
* Menu entries carry their state in the icon's **class** (`fa-check` / `fa-square`,
  `fa-dot-circle` / `fa-circle`); the `name` attribute stays `icon-square` even
  when ticked, so reading `name` reports everything as unticked. That cost a
  debugging cycle and is now documented in the spec.

## Automation limit found on dev (not filed)

**A context-menu submenu opens only once per page.** The first pointer move across
a group (`Statistics`, `Histograms`) opens its submenu; every later attempt in the
same page leaves it shut — verified with hover, with a click on the group, with the
pointer parked inside the popup first, and after dismissing the menu by clicking
away. Escape does not close these menus at all; only a click elsewhere does.

Consequences for this spec:

* exactly one submenu operation is automated (`Statistics > sum`);
* **Row Source** is driven through the Context Panel property grid instead of the
  All / Filtered / Selected radio items — same setting, reliable path;
* `Histograms > SEX`, removing a statistic again, and the date-column statistics
  stay manual.

## Not automated

* Histogram columns, removing a statistic, per-column visibility in the viewer's
  column manager, date-column statistics — `statistics.md` and `statistics-ui.md`.
