# Sunburst viewer — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI_v2.csv and demog.csv; add Sunburst for each | 9s | PASS | PASSED | `SPGI_v2.csv` not present on dev; used `SPGI.csv` (3624 rows) + `demog.csv` (5850 rows). `tv.addViewer('Sunburst')` adds the viewer cleanly in both cases. |
| 2 | Gear icon → Context Panel with viewer properties | 3s | PASS | PASSED | `sunburst.props.getProperties()` returned 8 properties across Data, Description, Misc, Color, Value. Asserted `hierarchyColumnNames`, `inheritFromGrid`, `includeNulls` present. |
| 3.1 | Table switching between SPGI and demog | 2s | AMBIGUOUS | FAILED | 2b used `closeAll` between file opens; live `sunburst.table = ...` rebind path was not exercised. Spec records as `test.skip(...)` inside softStep — outer stepErrors aggregator counts this as FAILED. |
| 3.2 | Select Columns dialog → choose 2-4 columns → OK | 3s | PASS | PASSED | `setOptions({hierarchyColumnNames: ['SEX','RACE']})` → readback via `props.get('hierarchyColumnNames')` matched. Select Columns dialog UI path and tooltip/label visibility not exercised. |
| 3.3 | demog: Inherit from grid with SEX coloring | 3s | PASS | PASSED | `setOptions({inheritFromGrid: true})` → readback `true`. Visual propagation of grid categorical coloring not verified. |
| 3.4 | SPGI: Core + R101 hierarchy; toggle Include nulls | 4s | PASS | PASSED | `setOptions({includeNulls: true})` then `false`, readback reflects each change. Null-segment rendering not verified visually. |
| 4 | Reset view (double-click empty space or context menu) | n/a | AMBIGUOUS | FAILED | Canvas-based; coordinate-precise double-click / context menu not synthesized. Recorded as `test.skip(...)` → FAILED. |
| 5 | Multi-selection (Click / Ctrl+Click / Ctrl+Shift+Click) | n/a | AMBIGUOUS | FAILED | Canvas-based selection; no sunburst selection API to assert result. Recorded as `test.skip(...)` → FAILED. |
| 6 | Select/filter on empty category (null segment) | n/a | AMBIGUOUS | FAILED | Canvas-based; depends on null-segment rendering and click. Recorded as `test.skip(...)` → FAILED. |
| 7 | Projects & layouts — save/restore | n/a | AMBIGUOUS | FAILED | Full layout round-trip with viewer-settings preservation is its own test surface. Recorded as `test.skip(...)` → FAILED. |
| 8 | Old layout compatibility — issue #2979 layout | n/a | SKIP | FAILED | Requires layout file from GitHub issue #2979 attachment, not available to the agent. Recorded as `test.skip(...)` → FAILED. |
| 9 | Collaborative filtering — internal + panel filters combine | n/a | AMBIGUOUS | FAILED | Not exercised in 2b. Recorded as `test.skip(...)` → FAILED. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 30s |
| grok-browser execution (scenario steps) | 12s |
| Execute via grok-browser (total) | 3m |
| Spec file generation | 1m |
| Spec script execution | 34s |
| **Total scenario run (with model)** | ~4m 34s |

## Summary

The sunburst viewer reproduces structurally on dev — Sunburst can be added to both SPGI and demog, and `hierarchyColumnNames`, `inheritFromGrid`, `includeNulls` set/readback works cleanly via `setOptions`/`props.get`. All interactive canvas steps (reset view, multi-selection, null-segment filtering) plus the layout round-trip, issue #2979 compatibility, and collaborative filtering were not exercised and are marked AMBIGUOUS/SKIP. Playwright run took 34s; `test.skip(...)` inside `softStep` is caught by the soft-step wrapper so every skipped step is reported as FAILED in the Playwright column even though the underlying behavior is "not exercised", not "broken". **Total scenario run (with model)**: ~4m 34s.

## Retrospective

### What worked well
- `addViewer('Sunburst')` followed by iterating `tv.viewers` to assert viewer presence is fast and stable.
- `sunburst.setOptions({...})` + `sunburst.props.get(...)` is a clean structural-verification path for every documented property (`hierarchyColumnNames`, `inheritFromGrid`, `includeNulls`).
- Splitting "open SPGI → open demog" into two `closeAll`-bounded evaluates kept runs deterministic across the two datasets.

### What did not work
- Scenario references `SPGI_v2.csv` which is absent on dev — only `SPGI.csv` and `SPGI_v2_infinity.csv` exist. Had to substitute `SPGI.csv`.
- Canvas-based interactions (double-click reset, Ctrl+Click multi-select, null-segment click) cannot be exercised without synthesizing coordinate-precise canvas mouse events and reading a selection API that Sunburst doesn't expose.
- `test.skip(...)` inside `softStep` throws a "Test is skipped" error that the softStep catcher treats as a step failure; AMBIGUOUS/SKIP steps thus appear as FAILED in Playwright column.

### Suggestions for the platform
- Expose a sunburst-viewer selection API (`currentSegment` / `selectedSegments`) so automated tests can verify click / Ctrl+Click / multi-select behavior without synthesizing canvas mouse events.

### Suggestions for the scenario
- `SPGI_v2.csv` is referenced but not present — use `SPGI.csv` or publish the `_v2` version to `System:DemoFiles`.
- Steps 4-9 require canvas-level interaction or external layout files; consider splitting visual/interactive steps from structural ones, or providing a helper that triggers them programmatically for regression coverage.
