# Delete — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Platform > Models | 30s | PASS | PASSED | Tree-click path failed — clicking `Predictive models` label inside Browse > Platform did not switch view (URL stayed `/`, view stayed `Home`). Fell back to JS API: `grok.shell.route('/models')` opened the Models browser (view.type=`models`, "0 / 0" empty state). |
| 2 | Find the model from the previous steps | 1m 45s | SKIP | FAILED | No model named `TestDemog` exists on dev for `agolovko`: `grok.dapi.models.list()` returned 0 total. Typed `TestDemog` into the search field — count stayed `0 / 0`. The prerequisite chain (train.md → browser.md / chemprop.md) was never completed on dev: `chemprop-run.md` recorded the chem-chemprop container wedged at `pending system stop`, and an MCP attempt to seed `TestDemog` via the Train Model UI was blocked by a Dart `NullError: method not found: 'aDm' on null` in `PredictiveModelingView.setModelNamePlaceholder` when SAVE was clicked. Playwright FAILED at `search.click()` actionability timeout 15s (the input is in the DOM but the Models view UI is in a transient state during the API call that returns 0). |
| 3 | Right-click it and select Delete | 0s | SKIP | FAILED | No model row to right-click. `test.skip()` inside `softStep` is caught and logged as a step failure. |
| 4 | In the confirmation dialog, click Delete | 0s | SKIP | FAILED | No dialog opens because step 3 was skipped. |
| 5 | Check that model has been deleted and is no longer present | 0s | SKIP | FAILED | Vacuous — there was nothing to delete. `grok.dapi.models.filter('name like "%TestDemog%"').list()` = 0 both before and after; cannot distinguish "deleted" from "never existed". |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 30s |
| grok-browser execution (scenario steps) | 45s |
| Execute via grok-browser (total) | 18m |
| Spec file generation | 4m |
| Spec script execution | 24s |
| **Total scenario run (with model)** | 22m 30s |

The bulk of "Execute via grok-browser (total)" was spent attempting the prerequisite — opening demog.csv, navigating ML > Models > Train Model, configuring Predict=SEX / Features=(11) All, expanding Misc/Description accordions, clicking SAVE twice, and confirming via console messages that a Dart NullError aborted the model-naming dialog.

## Summary

The delete UI itself could not be exercised on dev because the prerequisite predictive model from train.md/browser.md does not exist for `agolovko` on this server (`grok.dapi.models.list()` = 0). An MCP attempt to seed it failed with a Dart `NullError` in `PredictiveModelingView.setModelNamePlaceholder` that prevents the SAVE button from opening the model-naming dialog (xamgle/predictive_modeling_view.dart:580). Step 1 (navigation to /models) passed both in MCP and Playwright; steps 2–5 are SKIP in MCP and FAILED in Playwright (the spec asserts the prerequisite is met, then skips downstream steps). Total scenario run (with model): 22m 30s.

## Retrospective

### What worked well
- `grok.shell.route('/models')` reliably opens the Predictive models browser when tree navigation fails.
- The Models browser empty state (⊘ icon, "0 / 0" counter) is unambiguous — easy to detect "no models" both visually and via API.
- `grok.dapi.models.list()` and `grok.dapi.models.filter('name like ...').list()` give a clean machine-readable pre/post check that decouples verification from the DOM.

### What did not work
- Browse > Platform > Predictive models tree click is dead in the MCP run — clicking the `Predictive models` label inside the expanded Platform node neither selects nor opens the view; URL stays `/`. Same intermittent failure noted in `chemprop-run.md` step 4.1.
- Train Model SAVE on dev throws `NullError: method not found: 'aDm' on null` inside `PredictiveModelingView.setModelNamePlaceholder` (predictive_modeling_view.dart:580) — the model-naming dialog never opens, so no model can be persisted via the UI on dev. Console produces no user-facing balloon; the only signal is the Dart stack trace in the browser console.
- The Predictive Model view exposes no `input-host-Model-Engine` selector for demog/SEX target — engine selection is implicit and there is no way to override it from the form (chemprop saw the engine input only after Features was set, and only because Chemprop is the one engine for `Molecule`).
- The "Select columns" dialog for Features is canvas-rendered: typing into the Search field does not filter the list (the text field accepts input but the canvas keeps showing all rows), and clicking "All" inside an unfiltered dialog selects every column. There is no DOM-level handle for picking individual columns.
- `chem-chemprop` Docker container on dev is still wedged at `pending system stop` (re-confirmed via `chemprop-run.md` from the same week) — chemprop-based prerequisites are equally unrunnable on dev.

### Suggestions for the platform
- Fix the `setModelNamePlaceholder` NullError so SAVE always opens the naming dialog, or surface the failure to the user (balloon/toast) instead of swallowing it as a Dart console error.
- Make the Engine selector visible in the Train Model form even for non-chem datasets — users currently have no way to choose between EDA engines (XGBoost / Linear Regression / PLS / kNN) until after the model is auto-trained.
- Wire the "Select columns" dialog Search field to actually filter the canvas — typing should narrow the visible rows so "All" can be used to select a subset.
- Restart or replace the `chem-chemprop` container on dev — it has been stuck in `pending system stop` for at least the last week of test runs.
- Tree-click navigation to `Browse > Platform > Predictive models` should consistently open the Models view; today it silently no-ops in MCP/Playwright contexts and only routing works.

### Suggestions for the scenario
- Spell out the prerequisite explicitly: "Pre-condition: Train.md and/or Chemprop.md must have been completed in the same session, leaving at least one model owned by the current user." Today the line `Find the model from the previous steps` is interpreted as "any model" but the train chain may have failed silently.
- Name the model the scenario expects (e.g. `TestDemog`) so the search step has a definite target — same suggestion as the previous run's retrospective; still applies.
- Reconcile path wording across files: `delete.md` says `Browse > Platform > Models`, `browser.md` and `predictive-models.md` say `Browse > Platform > Predictive models`. The tree node is labeled "Predictive models".
- Add a final cleanup hint: "If the deletion succeeded, the search count goes from 1 / N to 0 / N-1 immediately; refresh the browser if the count does not update."
