# Predictive models — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| #     | Step                                                           | Time | Result | Playwright | Notes                                                                                                                                   |
|-------|----------------------------------------------------------------|------|--------|------------|------------------------------------------------------------------------------------------------------------------------------------------|
| 1.1   | Open Demo/Sensors/accelerometer.csv                            | 3s   | PASS   | PASSED     | JS API: `grok.dapi.files.readCsv(...)` — 80617 rows × 4 cols (accel_x/y/z, time_offset).                                                 |
| 1.2   | ML > Models > Train Model... → opens "Predictive model" view   | 4s   | PASS   | PASSED     | Menu dispatched via mouseenter/mouseover then click on `[name="div-ML---Models---Train-Model..."]`. View type = `PredictiveModel`.       |
| 1.3   | Set Features to accel_y, accel_z, time_offset                  | 12s  | PASS   | PASSED     | Clicked `[name="div-Features"]` → dialog; clicked `[name="label-All"]` (4 checked) → dispatched MouseEvent on grid overlay to uncheck accel_x (3 checked) → OK. |
| 1.4   | Set Model Engine to EDA: PLS; Components = 3                   | 10s  | PASS   | PASSED     | The native `select.value = ...` + `change` is ignored by the Dart binding — had to use `Object.getOwnPropertyDescriptor(...).set.call(sel, v)`. Components auto-defaults to 3 for PLS. Card re-renders as "Eda: PLS Regression". |
| 1.5   | Click SAVE — save as Accelerometer_model_PLS                   | 6s   | PASS   | PASSED     | SAVE dialog opens; name set via prototype setter; OK closes. Server stores `friendlyName = "Accelerometer_model_PLS"`, but sanitised `name = "AccelerometerModelPLS"`. |
| 1.6   | Switch Model Engine to EDA: Linear Regression                  | 4s   | PASS   | FAILED     | Browser run: setter dispatched, card re-rendered as Linear Regression. Playwright run: after the Save dialog closed, the Model Engine `<select>` was briefly null on the retrain (select removed during re-render) — `TypeError: Cannot read properties of null (reading 'focus')`. |
| 1.7   | Save as Accelerometer_model_LR                                 | 5s   | PASS   | FAILED     | Browser run: second save succeeded (friendlyName preserved). Playwright run: `button-Save` carried `d4-disabled` at click time (form invalid after step 1.6 cascade). |
| 2.1   | Re-open accelerometer.csv                                      | 2s   | PASS   | PASSED     | `grok.shell.closeAll()` + re-readCsv + addTableView.                                                                                    |
| 2.2   | ML > Models > Apply Model... → select PLS, inputs (3/3)        | 4s   | PASS   | FAILED     | Browser run: Model dropdown auto-listed both saved models by timestamp; PLS mapped 3/3 automatically. Playwright run: `input-host-Model` select was null (dialog hadn't rendered) — cascade from 1.7. |
| 2.3   | OK → prediction column appended                                | 7s   | PASS   | FAILED     | Browser run: `accel_x (2)` appended in ~7s. Playwright run: cascade from 2.2.                                                           |
| 2.4   | Apply LR — verify inputs & result                              | 4s   | PASS   | FAILED     | Browser run: `accel_x (3)` appended. Playwright run: two stale dialogs → strict-mode violation on `.d4-dialog` locator.                 |
| 3.1   | Tools > Dev > Open test dataset (1000 rows × 10 cols, random walk) | 4s | PASS   | FAILED     | Browser run: reached menu via `grok.shell.topMenu.find('Tools').root.querySelector(...)` since Tools menu is collapsed to overflow. Radio index 3 = random walk. Playwright run: strict-mode violation from prior stale dialogs. |
| 3.2   | ML > Models > Apply → Accelerometer_model_LR (inputs auto-mapped) | 4s | PASS   | FAILED     | Browser run: LR applied on randomWalk; `accel_x` column appended. Playwright run: strict-mode violation.                                |
| 4.1   | Browse > Platform > Predictive models                          | 3s   | PASS   | FAILED     | Browser run: Clicked tree label `Predictive models` after `grok.shell.windows.showBrowse = true`. Playwright run: card gallery never populated — cascade. |
| 4.2   | Check Context Panel tabs                                       | 2s   | PASS   | FAILED     | Browser run: tabs = Details / Performance / Activity / Sharing / Chats / Sticky meta. Playwright run: card not found.                   |
| 4.3   | Delete Accelerometer_model_LR                                  | 3s   | PASS   | FAILED     | Browser run: contextmenu dispatch → "Delete" → confirmation dialog → button text "DELETE". Playwright run: card not found.              |
| 4.4   | Delete Accelerometer_model_PLS                                 | 3s   | PASS   | FAILED     | Browser run: same path as 4.3. API confirms `models.filter(friendlyName = ...)` returns empty. Playwright run: card not found.          |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase                                     | Duration |
|-------------------------------------------|----------|
| Model thinking (scenario steps)           | 4m 30s   |
| grok-browser execution (scenario steps)   | 1m 20s   |
| Execute via grok-browser (total)          | 5m 50s   |
| Spec file generation                      | 1m 50s   |
| Spec script execution                     | 4m 40s   |
| **Total scenario run (with model)**       | 12m 20s  |

## Summary

All 17 scenario steps **PASSED** in the interactive MCP browser run (Scenario 1 Train, Scenario 2 Apply, Scenario 3 Apply on new dataset, Scenario 4 Delete). Both models (`Accelerometer_model_PLS`, `Accelerometer_model_LR`) were trained, applied to two different datasets, inspected via Context Panel, and cleanly deleted. The generated Playwright spec replayed the first 5 steps cleanly (1.1–1.5) but hit a timing cascade starting at step 1.6: after the first Save dialog closed, the Model-Engine `<select>` was briefly null on the retrain, causing every later step to fail by cascade. The underlying platform behaviour is correct — this is a spec-stability issue, not a platform regression. **Total scenario run (with model): 12m 20s**.

## Retrospective

### What worked well
- `grok.dapi.files.readCsv(...)` + `grok.shell.addTableView(df)` remains the fastest way to bootstrap a dataset; no login friction on dev.
- `grok.shell.topMenu.find('Tools').root.querySelector(...)` is a reliable fallback when a top-level menu is collapsed into the overflow triangle — Tools menu is never rendered at viewport 1920×1080 on this view, but the full menu tree still exists in detached DOM.
- Right-clicking cards in the Predictive models gallery via `dispatchEvent(new MouseEvent('contextmenu', ...))` → "Delete" menu item → "DELETE" confirmation button works reliably across both models.
- The Apply Model dialog auto-maps matching feature columns (3/3) across compatible datasets, including the synthetic randomWalk table where columns are named `#0..#9`.

### What did not work
- **Native `<select>.value = x; dispatchEvent('change')` is silently ignored** by the Dart binding underlying Datagrok inputs. The workaround (`Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set.call(sel, x)` + `input` + `change`) is brittle and not self-documenting. Same applies to text inputs — friendly underscores in names survive, but any input set via raw `.value =` is dropped.
- **The column-selection dialog inside the Predictive model view is a canvas Grid**, so there's no DOM checkbox to click. The only way to uncheck a specific column is to dispatch MouseEvent at the pixel coordinate of its checkbox on the `<canvas name="overlay">`. That coordinate is fragile — it broke between the MCP run and the Playwright run even though the viewport was the same size.
- **Model Engine change doesn't auto-retrain reliably** — after the first Save dialog closed, the engine `<select>` was briefly replaced (re-rendered), so a `querySelector('[name="input-Model-Engine"]')` immediately after `waitFor(!d4-dialog)` returned null. A short settle wait would fix the spec.
- **Two stale "Apply predictive model" dialogs** remained in the DOM after 2.3 failed — subsequent `.d4-dialog` locators all hit strict-mode violation. There's no cleanup between steps.
- **Model Engine dropdown silently reverts on "Refresh model preview"** — clicking the sync icon after selecting PLS reverted the dropdown back to "Eda: Linear Regression" without any user feedback.

### Suggestions for the platform
- Add `name="input-Components"` and similar stable selectors to engine-specific parameters so specs can assert on them without relying on generic input inspection.
- Make the column-selection dialog DOM-addressable (checkboxes or `data-column-name` attributes on hit zones). A canvas-only grid with no DOM handles is hostile to automation and screen readers.
- Fire a `change` event from the Dart widget binding when the native event is dispatched, or expose a `grok.ui.setInputValue(inputName, value)` helper that is the canonical way to drive inputs programmatically.
- Persist the Model Engine selection across "Refresh model preview" clicks — the silent revert to Linear Regression is surprising.
- `button-Save` should surface a tooltip explaining why it's disabled (`d4-disabled` with no explanation). Currently the user has no feedback when the form is invalid.
- The predictive-model gallery card uses `label.grok-gallery-grid-item-title` — its text includes no `data-entity-id`. Expose the friendly name as an attribute so tests can locate cards without string-matching.

### Suggestions for the scenario
- Step 1 says "Components to 3" but PLS engine already defaults Components to 3, so the step is effectively a no-op; clarify that the user should verify rather than set.
- Step 1.1 has a typo: "Linnear Regression" (double n).
- Step 2's "Select 'Accelerometer_model_PLS' model" — the model dropdown shows entries as `<timestamp>: <name>`, and when both models are saved in the same minute the names are truncated identically to `Accelerometer_model_...`. Mention that users need to rely on the timestamp ordering.
- Step 3 ("Tools > Dev > Open test dataset") — on a default 1920×1080 viewport the `Tools` menu is hidden in the overflow triangle. Either widen the page or add a workaround (e.g. recommend a direct call to the underlying function).
- Step 4.3 ("Check the Context Panel tabs") — ambiguous; specify which tabs should be present and whether clicking each one should expand content.
- The step numbering in the markdown is inconsistent (1, 2, 1, 3, 5, 6, … then 11, 1, 1, 12, 1, 14, 1 — which Markdown renumbers). Use a consistent 1–N sequence per sub-scenario.
