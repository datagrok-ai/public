# Radar viewer (Charts package) — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open earthquakes.csv; Add viewer > Radar | 3s | PASS | PASSED | `readCsv('System:DemoFiles/geo/earthquakes.csv')` + `addTableView` + semType wait + `addViewer('Radar')`. Viewers = [Grid, Radar]; rowCount=2426. No errors. |
| 2 | Open demog.csv; Add viewer > Radar | 3s | PASS | PASSED | `closeAll` + `readCsv('System:DemoFiles/demog.csv')` + `addTableView` + semType wait + `addViewer('Radar')`. Viewers = [Grid, Radar]; rowCount=5850. No errors. |
| 3 | Gear → Context Panel → check all properties (table, selection, Values, Style) | 1s | PASS | PASSED | Enumerated 21 properties via `radarViewer.props.getProperties()` across Data / Description / Misc / Selection / Color / Style / Value / Legend. Spec also toggled `backgroundMinColor` via `props.set` and confirmed getter echoed the new value. Live render interactivity not observed. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 45s |
| grok-browser execution (scenario steps) | 8s |
| Execute via grok-browser (total) | 1m 30s |
| Spec file generation | 35s |
| Spec script execution | 33s |
| **Total scenario run (with model)** | 2m 38s |

## Summary

The Radar viewer reproduced cleanly on dev for both earthquakes.csv (2426 rows) and demog.csv (5850 rows). All 21 Radar properties were enumerated across the required categories (Data, Selection, Value, Style, Legend), confirming the scenario's structural coverage items. The generated Playwright spec (self-contained, standalone login) passed all three `softStep()` sections in 33s. **Total scenario run (with model)**: 2m 38s.

## Retrospective

### What worked well
- `grok.dapi.files.readCsv` + `grok.shell.addTableView` + `tv.addViewer('Radar')` is a fast, deterministic path that bypasses the ribbon/gallery entirely.
- `radarViewer.props.getProperties()` + `props.set/get` made category verification and a Style-color roundtrip trivial — no gear-icon selector needed.
- Standalone Playwright spec (fresh context, no CDP attach) with the login snippet from SKILL.md worked first try once the `[name="Browse"]` wait was bumped to 120s.

### What did not work
- Scenario step 1 says just "earthquakes.csv" — actual dev path is `System:DemoFiles/geo/earthquakes.csv` (not top-level). Had to adjust the path at the MCP-run step.
- Live visual re-render on each property change was not observed (structural verification only); this is an inherent limitation of the JS-API-only path.

### Suggestions for the platform
- Consider adding `earthquakes.csv` as a top-level `System:DemoFiles/earthquakes.csv` alias or symlink so common references resolve without the `geo/` prefix.

### Suggestions for the scenario
- The path for earthquakes.csv on dev is `geo/earthquakes.csv` — the scenario should use a full `System:` path (e.g. `System:DemoFiles/geo/earthquakes.csv`) or move the file to the top `System:DemoFiles/` folder. The `datasets` JSON block should also list the earthquakes file.

---
{
  "order": 28,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/geo/earthquakes.csv"]
}
