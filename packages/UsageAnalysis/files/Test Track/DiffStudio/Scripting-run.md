# Verify Scripting and Model Interaction in Diff Studio — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1a | Open Diff Studio; turn on Edit toggle → equations editor opens | PASS | PASSED | Clicked `.ui-input-switch`; equations editor opened showing full Bioreactor IVP definition |
| 1b | Click </> icon on ribbon → JS script view opens | PASS | PASSED | Clicked uid of `</>` ribbon text; ScriptView opened with full JS source including //name:, //input:, //output: metadata |
| 2a | Execute the JavaScript script (Run button) | PASS | PASSED | Clicked `.grok-script-run-icon` play button; RichFunctionView opened with Inputs panel + Line chart + Grid |
| 2b | Move slider for Final at input; table and chart update in real time | PASS | PASSED | Changed Final range from 1000→500; Line chart x-axis shrunk to 500, curves updated. No Process mode input or Facet tab (only Multiaxis) as remarked |
| 3 | Add //tags: model to JS body; save script | PASS | PASSED | Used CodeMirror API to insert `//tags: model` after //meta.features line; Ctrl+S → "Script saved." toast appeared |
| 4 | Go to Apps > Run Model Hub; saved model appears in catalog | PASS | PASSED | Model Hub opened; Bioreactor visible in the 9-model catalog alongside GA-production, Pollution, PK-PD etc. |
| 5 | Interact with model in Model Hub; move slider updates table and chart | PASS | PASSED | Opened Bioreactor from Model Hub; changed Final from 1000→600; Line chart updated to t=600 in real time |

## Summary

All 5 steps (7 sub-steps) passed. The Edit toggle correctly opens the equations editor; the </> button opens the JavaScript script view. The script runs successfully producing the Bioreactor RichFunctionView. Adding `//tags: model` and saving registers the script in the Model Hub catalog, where it runs correctly with reactive slider updates.

## Retrospective

### What worked well
- Edit toggle (`.ui-input-switch` click) reliably toggled the equations editor on
- `</>` ribbon item click opened the JS script view in a ScriptView
- CodeMirror API (`cm.setValue()`) allowed reliable script editing
- Ctrl+S saved the script and showed a "Script saved." toast
- Model Hub correctly picked up the `//tags: model` script immediately after saving
- The saved model in Model Hub functioned identically to the DiffStudio version

### What did not work
- Navigating directly to `/apps/ModelCatalog` returned "Application not found" — correct URL is `/apps/Modelhub`
- The first approach to opening Model Hub (via `DG.Func.find({name: 'modelCatalog'}).apply()`) did not open a new view; direct URL navigation was needed

### Suggestions for the platform
- The `//tags: model` documentation could be more prominent in the script editor (e.g., a button/snippet for "Register as Model Hub model")
- The "Application not found" error banner persists even after a correct app is opened — this is a cosmetic UX issue

### Suggestions for the scenario
- The scenario says "Move the slider for the *Final at* input" but the input is labeled just "Final" in the UI (not "Final at")
- Step 4 should note the Model Hub URL is `/apps/Modelhub` (capital H, lowercase ub) for reproducibility
- Clarify that the model appears immediately in Model Hub after saving (no refresh needed)
