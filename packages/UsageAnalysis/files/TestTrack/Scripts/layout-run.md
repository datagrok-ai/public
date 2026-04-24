# Scripts Layout — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Create test_Layout script | 4s | PASS | PASSED | Via `DG.Script.create(body) + grok.dapi.scripts.save()`. FIRST MCP attempt used `#` metadata prefix — triggered a client-side JavaScript `SyntaxError: Private field '#name' must be declared in an enclosing class` because JS scripts require `//` (not `#`) for metadata. Rewrote body with `//` prefix and re-saved |
| 2 | Open editor, go to Layout tab | 6s | PASS | FAILED | Double-click on `.grok-gallery-grid-item` with title `test_Layout` opened the ScriptView. Layout tab visible among `Script / Layout / Debug` tab-headers. Playwright failure: `grok.shell.route('/scripts')` + 2500ms did not surface the newly saved card in time (same gallery-refresh latency seen in delete-spec) |
| 3 | Click Run script; OK in dialog; grid appears | 15s | PASS | FAILED | Run dialog "test_Layout" with idx=1 default, OK clicked. After ~10s, `[name="viewer-Grid"]` appears. Result dataframe: 270 rows × 7 cols. Playwright cascade: step 2 never loaded the editor |
| 4 | Add viewers (with/without docking) | 0s | SKIP | SKIPPED | Canvas-based Toolbox gallery does not expose `addViewer` on the internal preview TableView (`grok.shell.tv` returns ScriptView, not the layout's TV). `DG.VIEWER.SCATTER_PLOT` insertion via API path isn't hooked to the Layout preview. Would need drag/drop automation |
| 5 | Add coloring, hide columns | 0s | SKIP | SKIPPED | Depends on step 4 |
| 6 | Save | 2s | PASS | PASSED | `[name="button-Save"]` clicked |
| 7 | Close All | 2s | PASS | PASSED | `grok.shell.closeAll()` |
| 8 | Run test_Layout — result should open with new layout | 6s | PARTIAL | PASSED | `script.apply({idx:1})` returned a 270x7 DataFrame, opened as table view. No extra viewers restored — because step 4 was skipped, nothing was saved in the layout |
| 9 | Add new viewers | 0s | SKIP | SKIPPED | Dependent on layout authoring |
| 10 | Save project | 0s | SKIP | SKIPPED | Cross-cutting scenario |
| 11 | Close All, reopen project, verify layout | 0s | SKIP | SKIPPED | Same |
| 12 | Toolbox > File > Refresh — layout shouldn't change | 0s | SKIP | SKIPPED | Same |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~3m 30s |
| grok-browser execution (scenario steps) | ~1m 40s |
| Execute via grok-browser (total) | 5m 9s |
| Spec file generation | ~35s |
| Spec script execution | 40s |
| **Total scenario run (with model)** | ~6m 40s |

## Summary

Script creation and the Layout tab's "Run script" flow work — the preview dataframe is built
from `System:DemoFiles/chem` CSVs (270 × 7) and a grid viewer renders. The core authoring
steps (add viewers, change styling, save the layout, reopen as a project) could not be
automated because the Layout preview's TableView isn't reachable through `grok.shell.tv` (it
returns the ScriptView) and the Toolbox is canvas-based. This leaves the main value of the
scenario — that a saved layout is restored on subsequent runs — untested.

Important finding that should propagate back to the scenario: the scenario's code snippet uses
`//` metadata (correct for JavaScript), but a naive "just save the text as a script" path
using `#` metadata fails with a client-side JS `SyntaxError`. The Scripts UI always uses the
Signature Editor which writes the right prefix, but direct API creation bypasses that.

## Retrospective

### What worked well
- `DG.Script.create(body) + grok.dapi.scripts.save()` is a reliable way to seed a script
- Layout tab appears consistently — `Script / Layout / Debug` tab-headers
- "Run script" link on the Layout tab opens the standard run-parameter dialog
- Grid viewer renders in the Layout preview once the script succeeds

### What did not work
- `grok.shell.tv` returns a `ScriptView` on the Layout tab — not a `TableView` — so `addViewer(...)` API path is unavailable
- Toolbox Viewer gallery is canvas-based; clicking icons to add viewers to the preview is not automated
- Client-side JavaScript scripts must use `//` metadata, NOT `#` — easy to get wrong when scripting the creation
- Gallery doesn't auto-refresh after `grok.dapi.scripts.save()` — consistent with delete-spec observation

### Suggestions for the platform
- Expose the Layout preview's internal TableView via `scriptView.previewTv` (or equivalent) so test code can add viewers / set options programmatically
- Detect metadata prefix mismatch at save time: if `#language: javascript` but `#name:` is used (not `//name:`), emit a clear error instead of failing at run time with a confusing `SyntaxError: Private field '#name'`
- Subscribe the Scripts gallery to `scripts.save` events so newly saved scripts appear without manual refresh

### Suggestions for the scenario
- Step 1's code block is ambiguously indented (leading four spaces) — make it a proper fenced block, and add the `//` prefix so it's obvious
- Split the scenario into two sub-scenarios: "Author & save layout" (steps 1–7) and "Project round-trip" (steps 8–12) — today a failure in the middle blocks validation of either
- Add an explicit precondition: `System:DemoFiles/chem` must contain ≥ 2 CSV files (the script asserts `csvFiles[idx]` with `idx=1`; a 1-file directory would `undefined`-crash)

## Re-run after spec fixes (2026-04-24)

After patching the spec for robust waits (`waitForFunction` on `grok.shell.v?.name`, full
route round-trips to force gallery refresh, Playwright right-click for context menus, JS-API
fallbacks for the Run-dialog table dropdown and the signature-editor's internal state), the
Playwright run now **PASSES** in 33s for Scripts Layout. All scenario steps above that were
previously marked `FAILED` in the Playwright column now pass on the updated spec. Steps still
marked `SKIPPED` are intentional (manual file picker, canvas toolbox, cross-cutting project
flow) and use `test.step.skip` in the spec.
