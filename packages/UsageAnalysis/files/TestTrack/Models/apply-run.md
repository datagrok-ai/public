# Apply predictive model — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog.csv | 10s | PASS | SKIPPED | Loaded via `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')`. 5,850 rows, 11 cols. |
| 2 | Go to ML > Models > Apply Model... | 25s | PASS | SKIPPED | `[name="div-ML"]` clicked, submenu triggered by dispatching `mouseenter` on `[name="div-ML---Models"]`, then `[name="div-ML---Models---Apply-Model..."]` clicked. Dialog `[name="dialog-Apply-predictive-model"]` with title "Apply predictive model" opened. |
| 3 | Select the 'TestDemog' model | 85s | FAIL | SKIPPED | Model `<select>` in `[name="input-host-Model"]` never populated (0 options after 10s wait). `grok.dapi.models.list()` returned 0 models for the logged-in user on dev — TestDemog (trained on public in the sibling train.md scenario) does not exist on dev. Console logged recurring `No available adapters.` warnings. |
| 4 | Verify prediction column added to the table | 2s | SKIP | SKIPPED | Cannot apply a model that does not exist — column count stayed at 11 (unchanged from step 1); no prediction column added. |

**Time** = step 2b wall-clock (incl. thinking). **Result** = step 2b outcome. **Playwright** = step 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 40s |
| grok-browser execution (scenario steps) | 32s |
| Execute via grok-browser (total) | 2m 12s |
| Spec file generation | 30s |
| Spec script execution | 20s |
| **Total scenario run (with model)** | 3m 2s |

## Summary

Steps 1 and 2 passed — demog.csv opened and the "Apply predictive model" dialog opened correctly via the top menu. Step 3 failed because no "TestDemog" predictive model exists on dev.datagrok.ai for the current user (`grok.dapi.models.list()` returned an empty list and the dialog's Model dropdown stayed empty); the scenario's prerequisite (the sibling Models/train.md scenario) has not been run on this environment. Step 4 was skipped as a cascading consequence. The generated Playwright spec failed during the shared login helper (password input detached mid-focus on dev), so every softStep was SKIPPED before execution — this is a spec-login flake on dev, unrelated to the scenario logic. **Total scenario run (with model): 3m 2s.**

## Retrospective

### What worked well
- Opening demog.csv via `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` is fast and deterministic
- Top-menu navigation by `[name="div-ML"]` + `mouseenter` dispatch on `div-ML---Models` reliably expanded the submenu and revealed `div-ML---Models---Apply-Model...`
- The dialog exposes a clean `[name="dialog-Apply-predictive-model"]` root and inputs have `[name="input-host-<Caption>"]` wrappers — easy to target

### What did not work
- **TestDemog model is missing on dev** — the Apply scenario has an implicit prerequisite on Models/train.md, but models are not shared across environments, so running Apply standalone on dev always fails
- **Model dropdown gives no user feedback when empty** — the `<select>` silently shows nothing, with no "No models available" affordance or link to Train Model
- **Shared `loginToDatagrok` helper flaked on dev** — password input detached from the DOM mid-`focus()` call (`element was detached from the DOM, retrying`), timing out at 15s; this broke the whole spec before any scenario step could run

### Suggestions for the platform
- When the Model dropdown in the Apply dialog is empty, show placeholder text ("No models available — train one first") and a link to the Train Model dialog
- Surface the `No available adapters.` warning in the dialog UI (today it is only visible in the browser console) so users know why the dropdown is empty
- Consider listing models visible to the user's groups in the dropdown, not just the user's own models

### Suggestions for the scenario
- Step numbering has a typo (`1, 2, 1, 3` — should be `1, 2, 3, 4`)
- Add an explicit pre-condition line: "Run Models/train.md first on the same environment so that TestDemog exists"
- Clarify expected output in step 4 (e.g., "a new column named `SEX(2)` appears at the end of the table")
- Mention that the scenario will fail on any fresh environment where no model has been trained yet
