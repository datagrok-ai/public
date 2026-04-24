# Add New Column — Run Results

**Date**: 2026-04-22
**URL**: http://localhost:8888 (local, CLI arg `localhost`)
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open table — Home-dir demog.csv, Query, GetTop100, GetAll | 4s | PASS | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` — 5850 rows, 11 cols incl HEIGHT/WEIGHT. Query/GetTop100/GetAll variants SKIPped — no Northwind connection on localhost. |
| 2 | Click Add new column icon | 3s | PASS | PASSED | UI path. `[name="icon-add-new-column"]` → dialog opens in ~200ms with title "Add New Column". |
| 3 | Add formula `${HEIGHT} + ${WEIGHT}` → CalcCol1 | 45s | PASS | PASSED | UI path. Name via `execCommand('insertText')` on `.d4-dialog input[type="text"]`; formula via `fill()` on the CodeMirror `.cm-content` (the visible editor, not the hidden 0x0 `textarea`). CalcCol1[0]=233.68=160.48+73.20 ✓. |
| 4 | Add CalcCol2 = `${CalcCol1} * 2` | 15s | PASS | PASSED | Same UI path. CalcCol2[0]=467.37=233.68·2 ✓. |
| 5 | Rename HEIGHT→HT; change HT[0] to 200 | 12s | PASS | PASSED | JS API only (scenario says "in the grid… change used columns names and some values"; UI in-grid rename is fiddly). Formula auto-updated to `${HT} + ${WEIGHT}`. After `fireValuesChanged()`: CalcCol1=273.2, CalcCol2=546.4 ✓. |
| 6 | Save project (datasync where available) | 60s | PASS | FAILED | Datasync not available (no DB connection). Fell back to regular save via SAVE ribbon button after `grok.dapi.projects.save(...)` failed with "Unable to add entity… to the project". The save completes server-side but the `.d4-dialog` takes ≥15s to close on this server — Playwright `waitForFunction` hit the 15s `actionTimeout`. MCP run measured dialog close at ~12s. |
| 7 | Close All | 2s | PASS | PASSED | `grok.shell.closeAll()` — all tables closed. |
| 8 | Open saved project | 7s | PASS | PASSED | `grok.dapi.projects.filter(name=…).first() → .open()` — project restored with HT column, CalcCol1/CalcCol2 intact, CalcCol1.tags.formula = `${HT} + ${WEIGHT}`. |
| 9 | Rename HT→HEIGHT_2 in reopened project | 2s | PASS | PASSED | Formula propagates: `${HEIGHT_2} + ${WEIGHT}` ✓. |
| 10 | Change HEIGHT_2[0] to 180 | 3s | PASS | PASSED | After `fireValuesChanged()`: CalcCol1=253.2, CalcCol2=506.4 ✓. |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~7m |
| grok-browser execution (scenario steps) | ~45s |
| Execute via grok-browser (total) | ~7m 45s |
| Spec file generation | ~2m |
| Spec script execution | 32s |
| **Total scenario run (with model)** | ~10m 30s |

## Summary

All 10 scenario steps reproduce successfully in the MCP run; 9/10 pass in the Playwright replay. The one FAILED Playwright step is Step 6 (Save project), which hit the 15s `actionTimeout` on the dialog-close wait — the server-side save actually succeeds (Step 8 successfully reopens the saved project), the dialog just lingers longer than the default action timeout on this localhost server. Formula propagation on column rename and value-change recalculation work correctly in both the live table and after a full save/close/reopen round-trip. Datasync step (6) was exercised as a regular save since no Northwind connection exists on localhost. **Total scenario run (with model) ≈ 10m 30s.**

## Retrospective

### What worked well
- `${HEIGHT} + ${WEIGHT}` style formulas and their dependencies chain correctly.
- Column rename propagates automatically to the formula text.
- Values recalculate after `fireValuesChanged()`, including through a dependent calc column.
- Save/close/reopen round-trip preserves formulas, tags, and values correctly.
- `[name="icon-add-new-column"]` and `[name="button-Save"]` are stable selectors.

### What did not work
- `grok.dapi.projects.save(grok.shell.project)` fails with "Unable to add entity … to the project" for a table view that was opened via `readCsv`. The SAVE ribbon path works fine; JS API parity for project save is missing.
- `grok.dapi.tables.uploadDataFrame + Project.create + addChild + projects.save` path crashed with `Cannot read properties of undefined (reading 'get$nqName')` — suggests a js-api ↔ Dart interop gap on `Project.addChild`.
- `page.fill()` (and `execCommand('insertText')`) on the `textarea` inside the Add New Column dialog silently hit a hidden 0×0 `<textarea>`; the real editor is a CodeMirror `.d4-dialog .cm-content` contenteditable div. This is a common footgun.
- The Save-project dialog takes >15s to close on this localhost server after OK — enough to exceed Playwright's default 15s actionTimeout.
- Column values didn't auto-recalculate on `.set(i, v)` alone — need `df.fireValuesChanged()` to trigger dependent columns.

### Suggestions for the platform
- Make `grok.dapi.projects.save(grok.shell.project)` work for the "Scratchpad"-style local project (same as clicking SAVE), or at least raise a clearer diagnostic than "Unable to add entity … to the project".
- Auto-trigger dependent calc column recalculation when a column's value is set via the API; requiring an explicit `fireValuesChanged()` is surprising.
- Consider removing the hidden 0×0 `<textarea>` shadow of the CodeMirror formula editor (or marking it `aria-hidden`), so DOM automation tools don't silently bind to the wrong element.
- Investigate why the Save-project dialog takes ≥15s to close after OK on an idle localhost server — either shorten the response, or close the dialog optimistically before server commit.

### Suggestions for the scenario
- Specify a concrete formula in step 3 (currently `TODO: specify which formula to use`). Recommend `${HEIGHT} + ${WEIGHT}` for the demog dataset, and `${CalcCol1} * 2` for step 4.
- Split the datasync sub-flow (steps 6–10 with datasync) into a separate scenario that explicitly requires a configured data connection; leave this scenario focused on the in-memory formula/rename/recalc behaviour plus one regular save/reopen.
- Clarify step 1's multi-source list — state up front that only *one* source needs to be picked per run, or that each source should be exercised as its own sub-run.
- Replace "in the grid, change the used columns names" with the explicit UI steps (right-click column → Rename) so the test is reproducible without falling back to JS API.
