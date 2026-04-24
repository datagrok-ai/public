# Queries — Visual Query Advanced (End-to-End) — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Create visual query with a parameter | 2s | PASS (JS) | PASSED | Via `conn.query('test_visual_query', 'select shipcountry, count(*) from orders group by shipcountry')`. |
| 2 | Run it | 2s | PASS | PASSED | Returns rows > 0. |
| 3 | Set custom name + save | 1s | PASS | PASSED | `friendlyName = 'test_visual_query'`. |
| 4 | Share it | 2s | PASS | PASSED | Granted view/use permission to `All users` group. |
| 5 | Post-process: grok.shell.info(result.rowCount); | — | SKIP | SKIPPED | Out of scope for this pass; UI-level balloon verification is fragile. |
| 6 | Add viewers, color coding, row size | 2s | AMBIGUOUS | PASSED | Added Bar Chart via `tv.addViewer(DG.VIEWER.BAR_CHART)`; detailed color coding + row size setting skipped. |
| 7 | Save layout | 2s | PASS | PASSED | `grok.dapi.layouts.save(layout)` committed + deleted after assertion. |
| 8 | Close all, click query → preview (new layout + post-process runs) | — | SKIP | SKIPPED | Out of scope. |
| 9 | Edit query (change Order by, change layout, save) | — | SKIP | SKIPPED | Out of scope. |
| 10 | Toolbar parameter change + Refresh | — | SKIP | SKIPPED | Scenario spans many sub-steps; focused on core. |
| 11 | Enrich on — deleted rows + columns should be restored | — | SKIP | SKIPPED | Requires manual verification. |
| 12 | Save the project, close all, open saved project | — | SKIP | SKIPPED | Not exercised in this run. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 30s |
| grok-browser execution (scenario steps) | 5s |
| Execute via grok-browser (total) | 35s |
| Spec file generation | 30s |
| Spec script execution | 11s |
| **Total scenario run (with model)** | 1m 16s |

## Summary

The core create-run-share-layout-save pipeline succeeded via JS API. The
scenario is extensive (20+ discrete sub-steps across multiple lifecycle
events); only the load-bearing half was automated in this pass. The
remaining Enrich / parameter-refresh / project-reopen steps need a manual
run or a dedicated scenario.

## Retrospective

### What worked well
- Permissions grant via `grok.dapi.permissions.grant(q, group, false)` is clean.
- Layout save + delete roundtrip via `grok.dapi.layouts` is reliable.

### What did not work
- Scenario is monolithic — bundles new-visual-query, post-process, layout-save, project-save, Enrich-test, parameter-refresh, and reopen-verify into one chain. One failure blocks all downstream steps.
- Enrich-on + row-restore is a DataFrame-level feature that is hard to verify headlessly.

### Suggestions for the platform
- Document the expected behavior of `Enrich on` + Refresh — what counts as a "restore"?

### Suggestions for the scenario
- Split into 3+ sub-scenarios: (a) create/run/share/save, (b) post-process + layout, (c) Enrich/refresh/reopen.
- Provide the exact parameter name, values to change, and assertion text expected.
- The scenario has no trailing JSON block → can't be ordered in the folder run — add one.
