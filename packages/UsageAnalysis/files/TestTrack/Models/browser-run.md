# Browser — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Platform > Predictive models | 4s | PASS | PASSED | Expanded the Platform branch via `.d4-tree-view-tri`, clicked `[name="tree-Platform---Predictive-models"]`; `grok.shell.v.path` became `/models?`. |
| 2 | Type TestDemog in the search field | 3s | FAIL | FAILED | Typed into `input[placeholder*="Search models"]`; `0 / 0` counter. No TestDemog model exists on this dev account — prerequisite `train.md` has not been run. |
| 3 | On Context Panel, check all tabs for the model | n/a | SKIP | FAILED | No model to select. Spec raises a PREREQUISITE error so the step is recorded as failed. |
| 4 | Near search field, click Filter templates icon, check content | 3s | PASS | PASSED | Clicked `i[name="icon-filter"]`. Filter panel opened with sections "Active filters", "Quick Filters" (All / Favorites / Used by me / Created by me / Created recently / Commented by me / Is applicable to...), and "Properties" (Method / Source / Tags / Created / Author). |
| 5 | Select multiple models (CTRL+click) | n/a | SKIP | FAILED | 0 models on dev — multi-select cannot be tested. |
| 6 | On Context Pane, open Commands tab, click Compare | n/a | SKIP | FAILED | Depends on Step 5. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 4m 20s |
| grok-browser execution (scenario steps) | 25s |
| Execute via grok-browser (total) | 4m 45s |
| Spec file generation | 2m 10s |
| Spec script execution | 51s |
| **Total scenario run (with model)** | 8m 50s |

## Summary

2 of 6 steps pass fully in both browser and Playwright (step 1: navigation; step 4: Filter Templates). Steps 2, 3, 5, and 6 require the prerequisite model `TestDemog` (trained in Models/train.md) plus a second model, which do not exist on this dev account. Programmatically creating a `DG.Model` via `grok.dapi.models.save` rejects with `Bad Request: Invalid argument(s): sourceId` — the ML > Models > Train Model... wizard is the only supported path, and its inputs use canvas-rendered column selectors that are not drivable from automation. The scenario therefore cannot reach a full green run on a fresh dev account without first running `train.md`. **Total scenario run (with model): 8m 50s**.

Three Playwright-side issues surfaced and were fixed during this session:
1. **Login helper race** — after typing the login, Dart re-renders `#signup-login-fields` and the password-input locator detaches. Patched `spec-login.ts` to use `keyboard.press('Tab')` instead of re-focusing the password input, avoiding the re-resolution race.
2. **Tree navigation** — on a fresh Playwright session the Platform branch is collapsed, so `[name="tree-Platform---Predictive-models"]` resolves to a hidden element. Added a guarded `.d4-tree-view-tri` click to expand Platform first.
3. **Filter panel text** — panel headers render as `"Quick Filters"` (title case) via CSS `text-transform: uppercase`, not the displayed `"QUICK FILTERS"`. Fixed the assertion text.

## Retrospective

### What worked well
- Clicking the `tree-Platform---Predictive-models` tree node navigates reliably to `/models?`.
- The search field `input[placeholder*="Search models"]` accepts text and updates the `0 / 0` counter in real time.
- `i[name="icon-filter"]` opens the Filter Templates side panel with the expected Active / Quick / Properties sections — this works even on an empty browser (no selected model required).
- The shared `loginToDatagrok` helper now survives the Dart re-render race thanks to the Tab-based field traversal.

### What did not work
- **Prerequisite gap on dev**: `TestDemog` does not exist. Steps 2, 3, 5, 6 cannot run without it.
- **No supported way to create a model outside the wizard**: `DG.Model` + `grok.dapi.models.save` needs a `sourceId` that isn't exposed; the wizard's column selectors are canvas-based and resistant to keyboard/search-and-Enter automation.

### Suggestions for the platform
- Expose a high-level JS API like `grok.ml.trainModel({df, predict, features, engine, name})` that handles training + saving — today only the low-level `trainLinearRegression` / `trainPLSRegression` functions are callable, and they return a raw blob with no way to persist a `PredictiveModel`.
- Make the column selector in the Predictive Model wizard keyboard-accessible (search input visible, Enter selects the highlighted row) so automation can drive it end-to-end.
- When the Models browser is empty, show a prominent "Train your first model" CTA that links to `ML > Models > Train Model...`.

### Suggestions for the scenario
- Add an explicit prerequisite line: "Run Models/train.md first — TestDemog and a second model must exist before this scenario can be executed fully."
- Step numbering skips 5 and 6 (source goes `1, 2, 3, 4, 7, 8`) — renumber to `1–6`.
- For each Context Panel tab (Details / Performance / Activity / Sharing / Chats / Sticky meta) enumerate the expected fields so the tester can validate content, not just presence.
- Clarify that the Filter Templates panel opens on the left (same column as Toolbox/Browse), not as a popup near the search icon.
