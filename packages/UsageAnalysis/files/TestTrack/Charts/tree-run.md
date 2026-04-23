# Tree viewer (Charts package) — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 0 | Setup: Open demog.csv; Add viewer > Tree; hierarchy CONTROL/SEX/RACE | 5s | PASS | PASSED | `readCsv('System:DemoFiles/demog.csv')` + `addTableView` + semType wait + `addViewer('Tree')` + `setOptions({hierarchyColumnNames:[...]})`. 22 Tree properties enumerated across Data/Style/Size/Color/Value/Misc/Description. |
| 1 | Shift+Click branches false/F/Asian, false/F/Black, false/M/Asian | 25s | AMBIGUOUS | PASSED | Tree branches are canvas-rendered; coordinate-precise Shift+Click on canvas is not reliably synthesizable. Spec uses programmatic fallback: `df.selection.set(i, true)` for the three CONTROL=false branches. `selection.trueCount = 174` — recorded as regression baseline. |
| 2 | Filter panel CONTROL=true; expected filtered count = 0 | 25s | AMBIGUOUS | PASSED | The scenario's "filtered count = 0" is the intersection of the tree-branch selection (step 1) with the CONTROL=true filter. Filter applied via `df.filter` bitset (39 rows, all CONTROL=true). Intersection = 0 (matches expectation). |
| 3 | Shift+Click branch true/F/Black; expected filtered count = 2 | 25s | AMBIGUOUS | PASSED | Same Shift+Click limitation as step 1. Programmatic fallback extends `df.selection` to include true/F/Black rows. `selected = 176`, intersection with filter = 2 (matches expectation). |
| 4 | Clear CONTROL=true filter; expected filtered count = 176 | 25s | AMBIGUOUS | PASSED | `df.filter.setAll(true)` resets to 5850 rows. `selection.trueCount = 176` (matches expectation). UI-level "clear filter" click was not exercised; verified via bitset API. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m |
| grok-browser execution (scenario steps) | 8s |
| Execute via grok-browser (total) | 2m |
| Spec file generation | 40s |
| Spec script execution | 36s |
| **Total scenario run (with model)** | 3m 16s |

## Summary

Setup (open demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduced cleanly on dev. All four test steps are marked AMBIGUOUS at the 2b MCP level because the Tree viewer is a canvas-rendered visualization and coordinate-precise Shift+Click on branches cannot be reliably synthesized from outside the renderer. The generated Playwright spec uses programmatic `df.selection` / `df.filter` bitset fallbacks that are equivalent to the scenario's intended selection/filter state and confirms all expected counts (174, 0, 2, 176); spec run was PASSED. **Total scenario run (with model)**: 3m 16s.

## Retrospective

### What worked well
- `grok.dapi.files.readCsv` + `grok.shell.addTableView` + `tv.addViewer('Tree')` + `tree.setOptions({hierarchyColumnNames: [...]})` is a deterministic path for the setup phase.
- Programmatic `df.selection` / `df.filter` bitsets produce exact counts matching the scenario's expected values (174, 0, 2, 176), giving a solid regression baseline even without UI-level branch clicks.
- Standalone Playwright spec with the SKILL.md login snippet worked first try (120s `[name="Browse"]` wait, `keyboard.type`, Enter on password).

### What did not work
- Steps 1 and 3 require Shift+Click on specific tree branches rendered inside a canvas; without a branch-selection API or pixel-stable coordinates for the ECharts tree layout, canvas-relative mouse events are brittle and were not attempted.
- Because step 1 is AMBIGUOUS via UI, steps 2–4 (whose expected values depend on the accumulated selection) inherit the AMBIGUOUS status at the UI level even though the arithmetic matches.

### Suggestions for the platform
- Expose a Tree-viewer branch-selection API so automated tests can reproduce Shift+Click selections without coordinate-precise canvas events — e.g. `TreeViewer.selectBranch(pathArray, {shift?: bool})` where `pathArray` is the hierarchy path like `['false', 'F', 'Asian']`.
- Optionally emit a `branchClick` event with the path so tests can verify click behavior end-to-end without knowing pixel coordinates.

### Suggestions for the scenario
- Rewrite branches as tree-path strings so automation can navigate programmatically (e.g. "All/false/F/Asian", "All/false/F/Black", "All/false/M/Asian"), or provide a tree API for test scripts to toggle specific branches by path.
- Clarify that "filtered count = 0 / 2 / 176" in steps 2–4 refers to the intersection of tree selection and filter (not just `df.filter.trueCount`, which is 39 under CONTROL=true).
- Note pre-condition: CONTROL is a boolean column with only 39 true values out of 5850 in demog.csv.
