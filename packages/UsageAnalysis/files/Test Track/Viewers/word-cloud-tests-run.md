# Word Cloud tests — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Setup: close all, open demog, add Word Cloud | 25s | PASS | PASSED | Icon `[name="icon-Word-cloud"]` click adds viewer; canvas empty until column set. |
| 2 | Column assignment: RACE | 8s | PASS | PASSED | UI click on column selector did not open picker (silent no-op); used JS API `setOptions({column})`. Caucasian/Black/Other/Asian render. |
| 3 | Column assignment: DIS_POP | 4s | PASS | PASSED | Renders Indigestion, RA, Psoriasis, UC, AS, PsA. |
| 4 | Column assignment: SITE | 6s | AMBIGUOUS | PASSED | SITE column does not exist in demog.csv; substituted SEVERITY. |
| 5 | Column assignment: SEX | 4s | PASS | PASSED | Only F and M appear. |
| 6 | Shape: circle/diamond/triangle-fwd/triangle/pentagon/star | 18s | PASS | PASSED | `<select>` editor available; JS API accepted all 6 shape values. |
| 7 | Font size range 10/60, 20/20, 12/48 | 6s | PASS | PASSED | All values applied; single-size collapses range visually. |
| 8 | Text rotation 0/0, -90/90, step 45, defaults | 6s | PASS | PASSED | All rotation values applied. |
| 9 | Grid size 2, 20, 8 | 5s | PASS | PASSED | Values applied. |
| 10 | Draw out of bound on/off | 3s | PASS | PASSED | Toggle works. |
| 11 | Font family serif/monospace/sans-serif | 4s | PASS | PASSED | All families applied. |
| 12 | Bold on/off | 3s | PASS | PASSED | Toggle works. |
| 13 | Tooltip on hover (word + count) | 9s | AMBIGUOUS | PASSED | Tooltip appears with row count (e.g. "5267 rows") and updates on different words, but does not include the word itself — scenario expected "word and its row count". |
| 14 | Hover updates tooltip on another word | 2s | PASS | — | Row count varies by word; verified 5267 rows vs 354 rows. |
| 15 | Click word selects rows | 7s | PASS | PASSED | Clicking "Caucasian" selects 5267 rows, "Other" selects 354. |
| 16 | Click another word updates selection | 2s | PASS | — | Selection replaces the previous, not additive. |
| 17 | Click empty space clears selection | 6s | FAIL | — | Selection remains (354 rows) after click on empty top-left and bottom-right of canvas. Word Cloud does not clear selection on empty-space click. |
| 18 | Filter SEX=M updates word cloud | 8s | PASS | PASSED | Filtered rows: 2607; RACE counts update (Caucasian 2444, Other 75, Black 53, Asian 35). |
| 19 | Remove filter restores proportions | 3s | PASS | — | Filter cleared, 5850 rows. |
| 20 | Resize viewer by dragging border | 5s | SKIP | — | `dockManager.resize` not a public method; programmatic border drag not feasible in MCP. |
| 21 | Alt+F expand to full screen | 4s | PASS | PASSED | Viewer width 461 → 1920, height 913 → 993. |
| 22 | Alt+F returns viewer to normal | 3s | PASS | — | Restored to 461×913. |
| 23 | Right-click opens context menu | 5s | PASS | PASSED | Menu has Full Screen, Save as PNG, Properties..., shape/font submenus, etc. |
| 24 | Context menu has standard viewer options | 1s | PASS | — | Properties..., Save as PNG, Full Screen present. |
| 25 | Click Properties opens property panel | 3s | PASS | — | Property grid visible with Data/Misc sections. |
| 26 | Error state: column with > 500 unique values | 5s | PASS | PASSED | Viewer shows "The Word cloud viewer requires categorical column with 500 or fewer unique categories" when column=USUBJID. |
| 27 | Error state clears after setting RACE | 2s | PASS | — | Error gone, cloud renders. |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 4m 0s |
| grok-browser execution (scenario steps) | 2m 30s |
| Execute via grok-browser (total) | 6m 30s |
| Spec file generation | 1m 15s |
| Spec script execution | 17s |
| **Total scenario run (with model)** | 8m 2s |

## Summary

Word Cloud renders correctly and accepts all 7 documented properties (column, shape, font size range, rotation, grid size, drawOutOfBound, font family, bold) via `setOptions`. The viewer correctly selects rows on word click, responds to filter changes, supports Alt+F fullscreen toggle, exposes standard context menu entries, and shows a helpful error message for columns with > 500 unique categories. Two defects stand out: empty-space click does not clear selection, and the tooltip shows only the row count, not the word label. Total scenario run (with model): **8m 2s**.

## Retrospective

### What worked well
- `setOptions(...)` on the viewer is a fast, reliable way to drive every documented look property.
- Shape `<select>` editor is exposed cleanly via `select.property-grid-item-editor-spinner` — the property grid respects string-enum types.
- The 500-category guardrail produces a clear, actionable in-viewer error message.
- Alt+F fullscreen toggle is symmetric and works reliably.
- Filter interaction updates the viewer in the same frame as the filter change.

### What did not work
- **Clicking the column selector `[name="div-column-combobox-column"]` in the property panel does not open a picker popup.** The `.d4-column-selector` root receives the click but no dropdown/list appears — silent no-op. Had to fall back to `viewer.setOptions({column})`. Likely the Charts package Word Cloud doesn't wire up the standard column combobox popup.
- **Empty-space click on the Word Cloud canvas does not clear row selection.** Scenario step explicitly requires it; observed behavior keeps the previously selected word's rows selected.
- **Tooltip shows only `N rows`, not the word/category label.** Scenario expects "word and its row count"; current output omits the word.
- **Programmatic border-drag resize is not practical via MCP** — relied on Alt+F instead; step 1 of "Viewer resize" skipped.

### Suggestions for the platform
- Fix column-picker popup on Word Cloud property panel so UI-first automation and real users can change the column without the Properties dialog.
- Word Cloud tooltip should include the category label, e.g. `Caucasian\n5267 rows`. Currently only row count is shown.
- Empty-space click on the Word Cloud canvas should clear the dataframe selection (matches Bar Chart / Tree Map behavior).
- Consider exposing a public `viewer.resize(w, h)` method so automation and layout tools can programmatically resize without manipulating dock internals.

### Suggestions for the scenario
- The `demog.csv` dataset has no `SITE` column — either substitute an existing categorical column (e.g. `SEVERITY` or `DIS_POP`) or document the missing-column case explicitly.
- "Click empty space — selection is cleared" is currently a failing assertion; if intended behavior is undefined, rewrite step to acknowledge that.
- Tooltip step should clarify the exact expected content (word + count vs. count only); current wording is aspirational.
- "Error state (too many categories)" step is vague about the cardinality column — suggesting `USUBJID` (5,850 unique IDs in demog) makes the step deterministic.
- "Resize the viewer by dragging its border" is hard to automate; add an alternative wording like "via `dockManager` or Alt+F" for scripted runs.
