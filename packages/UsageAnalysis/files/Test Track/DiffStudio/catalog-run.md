# Catalog (DiffStudio) — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open DiffStudio; load PK-PD from Library | PASS | PASSED | Called `runDiffStudio()` to get full app with complete ribbon; clicked open model combo → Library → PK-PD; model loaded (5 cols: t, depot, centr, peri, eff; 1210 rows; 4 charts) |
| 2 | Click Save to Model Hub icon | PASS | PASSED | Found `.diff-studio-ribbon-save-to-model-catalog-icon` (fa-layer-plus) icon at x=715; clicked it; script saved to `grok.dapi.scripts` with tag "model"; verified via API: script "PKPD_1" saved with id ending in "d24e-6bc74a67c9db" |
| 3 | Go to Apps > Model Hub | PASS | PASSED | Navigated to `/apps/Modelhub`; Model Hub launched showing 10 models: GA-production, Template123, PK-PD ×2, Pollution, Bioreactor, PK-PD, Ball flight, Acid Production |
| 4 | Click Refresh icon | PASS | PASSED | Clicked `.fa-sync` icon; catalog refreshed (10/10 items, count unchanged as model was already present) |
| 5 | Run PK-PD model from catalog | PASS | PASSED | Double-clicked first PK-PD tile; DiffStudio launched with PK-PD model (5 cols, 1210 rows; Multiaxis + Facet tabs; 4 line charts showing cyclic dose patterns) |
| 6 | Modify inputs; verify results update | PASS | PASSED | Changed count from 10 to 5; rows updated 1210→605 in real-time; charts x-axis shortened from ~120h to ~60h; 5 dose cycles visible in all charts |

## Summary

All 6 steps passed. The full DiffStudio workflow — load from Library, save to Model Hub, run from catalog, modify inputs — works end-to-end. The key requirement was using `runDiffStudio()` (not `demoSimPKPD()`) to get the complete ribbon including the Save to Model Hub icon. Once the script was saved with tag "model", it appeared immediately in the Model Hub catalog.

## Retrospective

### What worked well
- `runDiffStudio()` opens the full app with complete ribbon including the `fa-layer-plus` "Save to Model Hub" icon
- Library menu navigation via the open model combo popup was reliable once properly triggered
- `grok.dapi.scripts` API confirmed the save succeeded (verifiable programmatically)
- Model Hub double-click to run PK-PD worked immediately
- Input modification (count 10→5) updated charts in real-time without needing a Run button

### What did not work
- `demoSimPKPD()` / `pkPdNew()` create demo views with only partial ribbon (no Save to Model Hub icon)
- The `runDiffStudio()` call opened with a "Template" model by default — needed to navigate to PK-PD via the Library menu
- No visual toast/notification for "Saved to Model Hub" was visible in the screenshot (DOM notification element not present for brief notifications); had to verify via `grok.dapi.scripts` API
- Refresh click (Step 4) returned same 10 items — the newly saved PKPD_1 script was already counted before the refresh because the list is loaded fresh on app open

### Suggestions for the platform
- `demoSimPKPD()` and other demo functions should launch a full DiffStudio app (with Save to Model Hub icon), not a limited demo view
- The "Saved to Model Hub" toast notification should persist longer (or use a persistent notification) so automated tests can capture it
- Model Hub should show a count badge after refresh if new models were added

### Suggestions for the scenario
- Step 1 should clarify whether to use `runDiffStudio()` (full app) or `demoSimPKPD()` (demo); only the full app has the Save to Model Hub button
- Step 4 (Refresh) could verify that the newly saved model appears in the catalog (e.g. by searching for the model name)
- Step 6 should specify which inputs to modify (e.g. "change count from 10 to 5") for reproducibility
