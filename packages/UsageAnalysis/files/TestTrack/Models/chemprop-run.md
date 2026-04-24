# Chemprop model — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1.1 | Open smiles.csv | 15s | PASS | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv')` → 1,000 rows × 20 cols; `canonical_smiles` detected as `Molecule`, `RingCount` as `int`. |
| 1.2 | ML > Models > Train Model | 10s | PASS | PASSED | Click `[name="div-ML"]` + `mouseenter` on `div-ML---Models` + click `div-ML---Models---Train-Model...` opens a full-page `PredictiveModel` **view** (not a dialog). |
| 1.3 | Set Predict to "Ring Count" | 15s | AMBIGUOUS | PASSED | Scenario says `Ring Count`; actual column is `RingCount`. Opened `.d4-column-selector` via `mousedown`, typed `RingCount` into the focused `.d4-column-selector-backdrop`. Result: `[name="input-host-Predict"] .d4-column-selector-column` textContent = `RingCount`. |
| 1.4 | Set Features to `canonical_smiles` — expect progress bar | 90s | FAIL | FAILED | Features input opens a canvas-rendered `[name="dialog-Select-columns..."]`. Typed `canonical_smiles` into search; clicking the row checkbox via `document.elementFromPoint(x, y).dispatchEvent(new PointerEvent(...))` at `(canvas.right-20, canvas.top+34)` with `buttons: 1` registered `1 checked` in MCP. OK → Chemprop engine auto-selected (`[name="input-host-Model-Engine"] select.value = 'Chem: Chemprop'`). **No progress bar rendered at any point** — the scenario's "a model starts training, check the progress bar" expectation fails: training only starts after clicking TRAIN/SAVE and filling the naming dialog. Playwright FAILED because the canvas-positioned PointerEvent hit a different element in the fresh 1920×1080 context (`received: 0`) — the canvas click is session-fragile even with PointerEvent. |
| 1.5 | After training, check interactive dashboard | 5s | FAIL | FAILED | The first TRAIN click opens the model-naming dialog immediately — no training, no dashboard. The current UI conflates TRAIN and SAVE into one action; the scenario's "After the training is complete, check the interactive dashboard" never occurs before naming. Playwright cascade: `input-host-Activation select` is null because 1.4 never activated the Chemprop engine. |
| 1.6 | Change Activation/Split_type/Epochs, click TRAIN | 5s | PASS | FAILED | Set Activation=`LeakyReLU`, Split-type=`scaffold_balanced`, Epochs=`5` via `<select>.value=` + `change` event and `<input>.value=` + `input`/`change` events; clicked `[name="button-Save"]`. Naming dialog opens. Playwright cascade: null selects from 1.4. |
| 1.7 | Save model as `test_chemprop` | 75s | FAIL | FAILED | Typed `test_chemprop` into `[name="input-Name"]`, clicked `[name="button-OK"]`. Polled `grok.dapi.models.filter('name = "test_chemprop"').list()` for ~75s — stayed at 0. `grok.dapi.models.list()` returned 0 total (clean user), so no model ever got persisted. Root cause: `chem-chemprop` container is stuck in `pending system stop`, so the training request has no running engine to answer. No balloon, no error shown to the user. |
| 1.8 | Change Metric to `auc` (`roc`), click TRAIN — expect balloon | 5s | FAIL | FAILED | Metric `<select>` options: `mse, mae, rmse, bounded-mse, bounded-mae, bounded-rmse, r2, binary-mcc, multiclass-mcc, roc, prc, accuracy, f1`. Scenario's `auc` does not literally exist; closest is `roc` (ROC-AUC, classification-only). Set to `roc`, clicked `[name="button-Save"]` — **no balloon**, the naming dialog opened again. The view does not validate metric against `dataset_type=regression`. |
| 2.1 | Close All and open smiles_only.csv | 15s | PASS | PASSED | `grok.shell.closeAll()` then `readCsv('System:DemoFiles/chem/smiles_only.csv')` → 1,000 rows × 1 col (`canonical_smiles`, Molecule). |
| 2.2 | ML > Models > Apply Model... — select test_chemprop | 45s | FAIL | FAILED | `[name="dialog-Apply-predictive-model"]` opened. Polled 30 × 500ms — `[name="input-host-Model"] select` stayed at 0 options because `grok.dapi.models.list()` returns 0 for this user on dev. Cancelled. |
| 2.3 | Verify prediction column "Ring Count (2)" added | n/a | SKIP | FAILED | Nothing was applied. |
| 3.1 | Browse > Platform > Dockers — locate chem-chemprop | 30s | PASS | FAILED | `showBrowse = true` + click `[name="Browse"]` + click `[name="tree-Platform"]` + click `[name="tree-Platform---Dockers"]` → `/dockers?`. Search `chem-chemprop` found 1 DOM match. Playwright FAILED: `[name="tree-Platform---Dockers"]` resolved to a hidden node — Platform row click did not actually expand the tree under Playwright's event pipeline. |
| 3.2 | Right-click chem-chemprop → Stop, then Run | 75s | FAIL | FAILED | Context menu opens correctly with items `Run / Stop / Restart / Revalidate image / Download context / Credentials... / ID / Grok name / Markup / Copy / Add to favorites`. Clicked `Stop` → polled `grok.dapi.docker.dockerContainers.filter('name = "chem-chemprop"').list()[0].status` 10 × 3s — stayed at `pending system stop`. Reopened menu, clicked `Run` → polled 10 × 3s — still `pending system stop`. Container is stuck on dev; the state machine does not advance. |
| 4.1 | Browse > Platform > Predictive models | 5s | PASS | FAILED | Clicked `[name="tree-Platform---Predictive-models"]` → `/models?` view opened. Playwright cascade from 3.1's tree-expansion failure. |
| 4.2 | Search test_chemprop, check Context Panel tabs | 5s | FAIL | FAILED | `grok.dapi.models.list()` returned 0. Nothing to click; no Context Panel tabs to verify. |
| 4.3 | Share the model | n/a | SKIP | PASSED | Spec step is a no-op when the Sharing tab isn't visible — passed vacuously; not a real verification. |
| 4.4 | Delete the model | n/a | SKIP | FAILED | No model exists; `grok.dapi.models.delete(...)` threw `no test_chemprop to delete`. |

**Time** = step 2b wall-clock (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 7m |
| grok-browser execution (scenario steps) | 6m |
| Execute via grok-browser (total) | 13m |
| Spec file generation | 2m |
| Spec script execution | 3m 45s |
| **Total scenario run (with model)** | 18m 45s |

## Summary

Only 3 of 16 sub-steps pass fully in both browser and Playwright (1.1 open `smiles.csv`, 1.2 open Train view, 2.1 open `smiles_only.csv`); 1.3 is AMBIGUOUS (column is `RingCount`, scenario says `Ring Count`). The scenario is blocked by the same two independent dev failures the prior run logged: (1) the `chem-chemprop` Docker container is wedged at `pending system stop` and does not transition on either Stop or Run — neither via the right-click menu nor over 30s of polling; (2) without the running container, Chemprop training silently never produces a model artifact, so `test_chemprop` is never persisted (`grok.dapi.models.list()` returns 0 for the current user), and everything downstream (Apply, Browse, Share, Delete) cascades. Two product-level issues re-surface independently of the container: (a) no progress bar is rendered during training (step 1.4 expectation), (b) the Predictive Model view does not validate `metric=roc` against `dataset_type=regression` — the naming dialog opens instead of a balloon (step 1.8 expectation). The MCP run found that the canvas-rendered column checkbox responds to `PointerEvent` dispatched via `document.elementFromPoint`, but transcribing the same coordinates into a fresh Playwright context still missed (0 checked), so the Playwright Train section cascades from 1.4. **Total scenario run (with model): 18m 45s.**

## Retrospective

### What worked well
- `grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv')` and `.../smiles_only.csv` open deterministically; semantic types (`Molecule`, `molregno`) are detected before the grid-canvas guard resolves.
- Top-menu navigation `[name="div-ML"]` → `mouseenter` on `[name="div-ML---Models"]` → click `[name="div-ML---Models---Train-Model..."]` is reliable in both MCP and Playwright.
- The Predictive Model view exposes every hyperparameter as a named input host (`input-host-Activation`, `input-host-Split-type`, `input-host-Epochs`, `input-host-Metric`, `input-host-Dataset-type`, …). `<select>` elements accept `.value = ...` + `change` event, so hyperparameter tuning is automation-friendly without dialogs.
- The `PredictiveModel` Chemprop engine auto-selects the moment a `Molecule` column lands in `Features` (observed via `[name="input-host-Model-Engine"] select.value === 'Chem: Chemprop'`).
- The Dockers view (`/dockers?`) and its right-click menu (`Run / Stop / Restart / Revalidate image / …`) are driven entirely by DOM events — no canvas interaction needed for the menu itself.
- `PointerEvent` dispatched on `document.elementFromPoint(x, y)` with `buttons: 1` toggles a canvas-grid checkbox in MCP (where `MouseEvent` alone silently drops). Useful technique for canvas-based widgets, though the pixel target remains brittle across viewport changes.

### What did not work
- **`chem-chemprop` container stuck on dev**: status `pending system stop` for the full run, across both Stop-then-Run cycles and >30s polling. Training never runs.
- **Canvas-rendered Select columns dialog is session-fragile**: the PointerEvent pixel technique that works against the live MCP session still missed in the fresh `test({page})` Playwright context (0 checked). The canvas row-1 checkbox at `(canvas.right-20, canvas.top+34)` depends on the dialog's exact position, which shifts between sessions.
- **Silent training failure, no user feedback**: the naming dialog accepts `test_chemprop`, closes, and no balloon/toast/notification appears. The server-side `FileSystemException` on the missing `.bin` file never surfaces in the UI.
- **No progress bar during training**: step 1.4 expects visible training progress. The Predictive Model view shows nothing — no progress bar, no epoch text, no activity indicator.
- **No metric/dataset_type validation**: setting `metric=roc` against the default `dataset_type=regression` does not produce a balloon; the naming dialog opens as if the config were valid (step 1.8 expectation fails).
- **Tree-Platform expand logic in fresh Playwright context**: clicking `[name="tree-Platform"]` in a freshly logged-in session does not reliably expand the node — subsequent `[name="tree-Platform---Dockers"]` / `tree-Platform---Predictive-models` stay hidden. The heuristic needs a dedicated expand signal rather than a row-level click.
- **Scenario conflates TRAIN and SAVE**: scenario text uses "click TRAIN" then "save model with name …" as separate steps. The current UI has a single button labelled visibly as `TRAIN` but carrying `[name="button-Save"]`; clicking it opens the save dialog, and training starts on OK. The expected two-step flow does not exist.

### Suggestions for the platform
- Pre-flight the TRAIN action: before submitting, verify (a) the engine's Docker container status is in a startable state and (b) the `metric` is valid for the chosen `dataset_type`. On failure, show an inline balloon with a specific fix ("Metric 'roc' requires dataset_type=classification" / "Chemprop container is stopped — start it in Platform > Dockers").
- Surface container-level failures (`FileSystemException`, container not started) as user-visible balloons on the Predictive Model view, not just in the browser console.
- Render a training progress indicator — progress bar, per-epoch status, or an explicit "waiting for engine" state — so users can distinguish running vs. stuck vs. failed.
- Expose a first-class JS API for end-to-end train+save (`grok.ml.trainModel({df, predict, features, engine, hyperparameters, name, ...})`) so smoke tests and automation don't need to drive the canvas-backed Select columns dialog.
- Replace the canvas-based column grid in the Select columns dialog with a DOM-backed virtualized list of checkbox rows, keyboard-navigable. Alternatively, make the existing `Search` input actually filter the canvas grid to matching rows (currently it's cosmetic — "All" still checks every column regardless of the filter).
- Auto-recover or force-unstick containers in `pending system stop` after a timeout; add a "Force restart" action to the Dockers context menu.
- Publish stable JS hooks on `.d4-tree-view-node` so `[name="tree-Platform---Dockers"]` can be directly navigated/expanded without relying on row-click toggling.

### Suggestions for the scenario
- Add a **Pre-conditions** section: "Requires `chem-chemprop` Docker container to be in a startable state. If not, start it via `Browse > Platform > Dockers` before step 1 (see section 3 for the recipe)."
- Fix the column reference: "Ring Count" → `RingCount` (match the demo CSV).
- Rewrite the "click TRAIN, then save model" phrasing to match the current single-button flow: "Click **TRAIN/SAVE**; in the naming dialog, enter `test_chemprop` and click OK. Training begins once the dialog is confirmed."
- Replace `auc` with `roc` in step 1.8 (the UI does not expose an `auc` option) and state the expected balloon text explicitly: "Balloon: 'Metric roc requires dataset_type=classification'".
- Renumber the repeated `5`s in section 1 (current file has `1,2,3,4,5,5,5,6`).
- Section 4: list the expected Context Panel tabs (Details / Performance / Activity / Sharing / Chats / Sticky meta) and the expected column delta for the Apply step ("smiles_only.csv goes 1 → 2 columns; new column `Ring Count (2)` of type `double`").
