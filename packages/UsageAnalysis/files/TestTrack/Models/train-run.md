# Train — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| #  | Step                                                              | Time | Result    | Playwright | Notes                                                                                                                                                                |
|----|-------------------------------------------------------------------|------|-----------|------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1  | Open demog.csv                                                    | 4s   | PASS      | PASSED     | `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` → 5850 rows × 11 cols.                                                                                       |
| 2  | ML > Models > Train Model... → opens Predictive model view        | 6s   | AMBIGUOUS | PASSED     | Scenario says "A dialog opens" but a `PredictiveModel` view opens (not a modal). Hover via dispatched mouseenter on `[name="div-ML---Models"]` then click submenu.   |
| 3a | Set Predict to SEX                                                | 12s  | PASS      | PASSED     | Mousedown didn't open the popup — only the chrome-devtools `click` (real CDP input event) did. Then `type_text "SEX"` + Enter. Caption updated to SEX.                |
| 3b | Set Features to WEIGHT and HEIGHT                                 | 18s  | PASS      | PASSED     | Click on Features opens `dialog-Select-columns...` whose grid is canvas-only. No DOM checkboxes — dispatched MouseEvent at canvas pixel coordinates (x=826, y=397/425). |
| 4  | Select Impute missing — k-NN dialog opens                         | 3s   | PASS      | PASSED     | Click checkbox → `dialog-k-NN-Imputation` opens. Inputs: Impute, Using, Distance, Neighbors, In-place, Keep-empty.                                                   |
| 5  | Click RUN — model trains and result renders                       | 25s  | PASS      | PASSED     | `[name="button-RUN"]`. After ~10s: Eda XGBoost trained. Dashboard: scatter plot Predicted SEX vs Actual, Confusion matrix (2919/324, 439/2168), Acc=0.870.            |
| 6  | Unselect Impute missing                                           | 1s   | PASS      | PASSED     | Click toggled checkbox to false.                                                                                                                                     |
| 7  | Select Ignore missing                                             | 4s   | PASS      | PASSED     | Auto-retrained: Eda XGBoost, Specificity 0.836, Sensitivity 0.904, Accuracy 0.875.                                                                                   |
| 8  | Select Predict Probability — verify result                        | 4s   | PASS      | PASSED     | Engine switched to Eda: Linear Regression. Positive class cutoff 0.50 appeared. ROC Curve + Predicted SEX-actual scatter rendered.                                   |
| 9  | Save model as TestDemog                                           | 7s   | PASS      | PASSED     | SAVE button opens dialog (empty title, has Name/Description/Tags). Focus + type "TestDemog" → OK. `grok.dapi.models.filter('friendlyName = "TestDemog"')` → 1.        |
| 10 | Repeat: train regression — WEIGHT by HEIGHT, save as TestDemog2   | 35s  | PASS      | PASSED     | Predict→WEIGHT, Features→only HEIGHT (uncheck WEIGHT after dialog reorders checked rows to top), re-tick Ignore missing (was cleared on input change). MSE=329.88, RMSE=18.16, R²=0.19. Saved as TestDemog2. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase                                     | Duration |
|-------------------------------------------|----------|
| Model thinking (scenario steps)           | 4m 20s   |
| grok-browser execution (scenario steps)   | 1m 50s   |
| Execute via grok-browser (total)          | 6m 10s   |
| Spec file generation                      | 1m 50s   |
| Spec script execution                     | 28s      |
| **Total scenario run (with model)**       | 9m 10s   |

All rows are full-phase wall-clock (incl. model thinking and retries), not just tool latency.

## Summary

All 10 scenario steps passed in the interactive MCP browser run. Both models — classification (TestDemog, Predict probability + Linear Regression) and regression (TestDemog2, Linear Regression on WEIGHT by HEIGHT) — trained successfully on dev.datagrok.ai. The generated Playwright spec replayed the full flow in 28s in its first run; an observed second-run flake on the regression step (save dialog from the classification flow occasionally remains in the DOM after OK click and intercepts subsequent clicks) is documented below as a stability note rather than patched, per the skill's no-second-editing-window rule. **Total scenario run (with model): 9m 10s**.

## Retrospective

### What worked well
- The Predictive Model **view** (not a dialog) is laid out with stable `[name="input-host-*"]` selectors for every input (Predict, Features, Ignore missing, Impute missing, Predict probability, Model Engine, etc.) — easy to introspect and assert on.
- Selecting "Impute missing" cleanly opens `dialog-k-NN-Imputation` with `[name="button-RUN"]` and `[name="button-CANCEL"]` — both have stable `name=` attributes.
- `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` + `grok.shell.addTableView(df)` is the fastest reliable bootstrap.
- The model-save flow integrates cleanly with `grok.dapi.models.filter('friendlyName = "..."').list()` — perfect for asserting persistence in tests.
- Real Chrome-DevTools `click` events (CDP input synthesis) reliably open the column-selector popup; synthetic `MouseEvent`/`PointerEvent` dispatched via `evaluate_script` are silently ignored by the underlying Dart widget.

### What did not work
- **Synthetic mousedown on the column-selector caption is dropped.** `dispatchEvent(new MouseEvent('mousedown', ...))` on `[name="div-Predict"]` does not open the column popup, even though `column_combo_box.md` says the popup opens on mousedown. The CDP `click` tool works because it goes through Chrome's input pipeline; raw JS dispatch does not. This is the same Dart change-listener gap that bites text inputs.
- **The Select-columns dialog is a fully canvas-rendered grid** — no DOM checkboxes — so unchecking a specific column requires dispatching MouseEvent at hard-coded pixel offsets on `canvas[name="overlay"]` (~y=257 + 28*rowIndex, x=826 for the checkbox column). The coordinates are fragile and the row order changes when the dialog reopens (checked rows move to the top), so the same pixel hit toggles a different column on the second visit.
- **Save dialog occasionally persists after OK click** — the `[name="dialog-"]` save dialog (empty `name` because `name=""` from an unset title) sometimes does not detach from the DOM after the OK click even though the model was saved on the server (verified via `grok.dapi.models.filter`). On the next step, the lingering dialog intercepts pointer events and blocks all further clicks. Observed on the second Playwright run; first run was clean.
- **Toggling Impute missing → triggers k-NN dialog every time, even on uncheck/recheck**, requiring an extra dismiss action that's not in the scenario script.
- **Switching the target column resets all "actions" checkboxes silently** (`for (BoolInput checkbox in actionsCheckboxes) { checkbox.value = false }` in the Dart view), so after switching to a regression target the user must re-check Ignore missing before training resumes — easy to miss.

### Suggestions for the platform
- Open the column-selector popup on a generic interaction event (e.g. `pointerdown` listener that also fires on Dart-internal synthetic dispatch), or expose a `grok.ui.openColumnSelector(input)` helper. Right now nothing short of real CDP input drives it.
- Make the `dialog-Select-columns...` grid DOM-addressable: add `data-column-name` attributes on hit zones (or render checkboxes as real `<input type="checkbox">` overlays). A canvas-only grid with no DOM handles is hostile to automation, screen readers, and keyboard navigation.
- Stabilise the row order in the Select-columns dialog — promoting checked rows to the top makes coordinate-based clicks stale and breaks any spec that doesn't re-introspect the dialog state between visits.
- Set a stable `name=` on the model-save dialog (e.g. `dialog-Save-Predictive-Model`) so specs don't have to scope by `.d4-dialog [name="input-host-Name"]`. The empty `name="dialog-"` is also a hygiene issue.
- Ensure the save dialog is reliably removed from the DOM after `button-OK` (await the underlying save before tearing down the modal, or tear down on the click handler then await, but don't leave it half-attached). Add a `dialog.close()` in a `finally` if the OK handler can throw.
- Consider _not_ resetting all action checkboxes when only the target column changes — at minimum, Ignore missing is a per-dataframe concern, not a per-target concern.

### Suggestions for the scenario
- "**A dialog opens**" in step 2 is wrong — Train Model opens a **view** (`Predictive model`), not a dialog. Scenario should say "A predictive-model view opens with Predict / Features / option checkboxes".
- Step 4's "**a dialog opens**" is correct, but should clarify the dialog title (`k-NN Imputation`) and the engine that gets used (it requires the EDA package's k-NN imputation, which depends on a non-trivial Docker container on some environments).
- Step 5 ("Click RUN — check the result") is ambiguous — what should the user verify? Suggest: "Verify the dashboard renders with a confusion matrix and accuracy ≥ 0.85".
- Step 8 ("Predict Probability — check the result") is similarly vague — specify that the engine should switch to a regression engine, ROC Curve should appear, and a "Positive class cutoff" input should be visible.
- Step 9 ("Save the model as TestDemog") doesn't mention that the dialog also has a Description and Tags field, nor that the OK button stays enabled even when the form is technically invalid (placeholder is auto-filled).
- The "Repeat training the model predicting WEIGHT by HEIGHT" line should explicitly enumerate the substeps, especially that the user must (a) deselect WEIGHT from Features (since the predict column can't also be a feature), and (b) re-enable Ignore missing because it's auto-cleared on target change.
- The numbered list uses `1.` repeatedly (Markdown re-numbers) — switch to a continuous 1–N sequence so step numbers in the rendered scenario match what tests/specs reference.
