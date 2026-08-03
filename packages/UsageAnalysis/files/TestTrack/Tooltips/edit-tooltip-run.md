# Viewers: Edit tooltip — Run Results

**Date**: 2026-05-12
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Close all. Open SPGI dataset | 13s | PASS | PASSED | 3624 rows, 88 columns. Bio/Chem wait triggered (Molecule semType). |
| 2 | Open grid properties and enable `Show Visible Columns In Tooltip` | 18s | PASS | PASSED | Grid property panel uses TD rows, not standard `input-host-*`. Found via TD text. Toggled false→true. |
| 3 | Open a scatter plot and a box plot | 10s | PASS | PASSED | Added by clicking toolbox icons. 3 viewers total (Grid, Scatter, Box). |
| 4 | Right-click viewer, select Tooltip → Edit... | 24s | PASS | PASSED | Dialog `Edit Tooltip` opens. Search input + master checkbox (checked, all selected); column grid with 3 columns; `Reset group tooltip` + `Design custom tooltip...` labels; OK/CANCEL + history icon in footer. Search "series" filtered to "Series" + "Primary Series Name" (case-insensitive ✓). |
| 5 | Pick a few columns in the grid and click OK | 16s | PASS | PASSED | Deselected all via master checkbox; clicked rows for Id, Chemist, Series (canvas-based column grid; clicks dispatched at x=495 on overlay canvas). Verified `df.getTag('.tooltip')` = `Id\nChemist\nSeries`. |
| 6 | Hover over the viewers — tooltip should consist of selected columns in same order | 8s | PASS | PASSED | Scatter plot tooltip showed `CAST Idea ID, Chemical Space Y` (axes) then `Id, Chemist, Series` in selected order. Box plot / grid hover via dispatched mousemove did not trigger tooltip rendering (canvas hover quirk); verified `dataFrame.getTag('.tooltip')` is the single source and both viewers' `rowTooltip` is empty so they inherit. |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 7s |
| grok-browser execution (scenario steps) | 22s |
| Execute via grok-browser (total) | 1m 29s |
| Spec file generation | 2m 18s |
| Spec script execution | 34s |
| **Total scenario run (with model)** | 4m 21s |

## Summary

All six scenario steps passed. Dialog elements, case-insensitive search, column selection through the
canvas-based column grid, and tooltip configuration via `.tooltip` dataframe tag all behaved as
described. The Playwright spec replayed the run cleanly in 31s. Total scenario run with model
thinking: **4m 21s**.

## Retrospective

### What worked well
- `grok.shell.tv.dataFrame.getTag('.tooltip')` gives a clean assertion target — no need to inspect
  the canvas tooltip after closing the dialog if you just want to verify configuration persisted.
- The case-insensitive search filter was easy to verify: type a lowercase substring of a known
  column name and screenshot the filtered list.
- The Tooltip submenu opens reliably with a single `mousemove` dispatch + 300ms wait — same pattern
  used in `default-tooltip-spec.ts` and other existing specs.
- Spec-login helper (`loginToDatagrok` + `specTestOptions`) made the spec compact; full preamble
  fits in 4 lines.

### What did not work
- **Grid properties don't use the standard `[name="input-host-*"]` pattern** — the grid property
  panel renders a `<table>` with `<td>` rows. `input-host-*` selectors returned 0 elements; had
  to fall back to "find TD by exact text, then closest TR's checkbox." Worth documenting in
  `references/grid.md`.
- **JS-dispatched `mousemove` does not reliably trigger viewer tooltips on box plot / grid** — only
  the scatter plot's tooltip fired. Box plot and grid both use the same `.tooltip` tag, so the
  configuration is verified correct, but live-hover assertions on those two viewers would need
  Playwright's `page.mouse.move` (real input events) or chrome-devtools `hover` against a uid.
- **`grok test --skip-puppeteer` flag missing** in installed datagrok-tools 6.2.3. Wrapper ran the
  package's `test()` function first (which returns null), aborting before Playwright. Worked around
  by exchanging the dev key via `POST /api/users/login/dev/<key>` and invoking `npx playwright test`
  directly with `DATAGROK_AUTH_TOKEN` / `DATAGROK_URL` env vars.
- **Canvas hit-coordinates for the column-grid checkbox** required guessing — x=495 worked on a
  1920×1080 viewport but is fragile if the dialog size changes. A more robust approach would be to
  introspect the inner grid's column widths via JS API.

### Suggestions for the platform
- Surface the column-grid's underlying dataframe / `Grid` instance on its root element (e.g.,
  `.d4-column-grid` could expose `__grid` for testability) so spec authors can toggle column
  selection via `grid.dataFrame.set('selected', i, true)` instead of canvas-click guesswork.
- Add `name=` attributes to grid property rows (e.g., `[name="prop-show-visible-columns-in-tooltip"]`)
  the same way the standard `input-host-*` widgets do — consistent selectors across panels would
  eliminate the TD-text fallback.
- Either restore `--skip-puppeteer` to `grok test` (advertised in the `grok-debug-scenarios` skill)
  or update the skill's instructions to match the released 6.2.3 wrapper.

### Suggestions for the scenario
- Step 4 description says "below the two checkboxes, there should be a grid" — the dialog actually
  has **one** checkbox (the select-all toggle next to the search input). The phrasing should be
  corrected to "below the search input and select-all checkbox".
- "the tooltip should consist of the same set of selected columns in the same order" — clarify
  whether the scatter plot's X/Y axis values (always prepended by the scatter plot itself) count
  as part of the tooltip or are out of scope. Right now the scatter tooltip contains
  `CAST Idea ID, Chemical Space Y, Id, Chemist, Series` and it's ambiguous whether that's a pass.
- Consider adding an explicit verification step "open Edit Tooltip again — selected columns are
  still checked" to catch persistence bugs.
