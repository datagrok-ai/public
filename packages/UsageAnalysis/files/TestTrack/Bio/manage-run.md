# Bio Manage Monomer Libraries — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open sample_Helm.csv | 12s | PASS | PASSED | Loaded `System:AppData/Bio/samples/HELM.csv` via `grok.dapi.files.readCsv` (5.6s incl. semType wait + Bio settle). 540 rows; HELM column → `Macromolecule`, Activity → no semType. Scenario file name "sample_Helm.csv" is a misnomer; actual dataset is `HELM.csv`. |
| 2 | Go to Bio > Manage > Monomer Libraries | 6s | PASS | PASSED | Clicked `[name="div-Bio"]`, dispatched `mouseenter`+`mouseover` on `[name="div-Bio---Manage"]`, clicked `[name="div-Bio---Manage---Monomer-Libraries"]`. Opens a **View** (`grok.shell.v.type === "view"`, name "Manage Monomer Libraries"), not a modal dialog. ~1.1s of MCP latency. |
| 3 | Check dialog output and checkboxes functionality | 18s | PASS | PASSED | Found 5 checkboxes — `HELMCoreLibrary.json`, `HELMCoreLibrary123456.json`, `NH2.json`, `polytool-lib.json`, `sample-lib-Aca-colored.json` — all checked by default. Toggling `HELMCoreLibrary123456.json` off shrinks the duplicate-symbols panel (9 → 3 cards; 303 → 129 child elements) and fires balloons: `Monomer library user settings saved` + `Monomer lib updated: PEPTIDE 336RNA 383`. Re-toggling restores cards (3 → 9) and child count (129 → 303). No stray error balloons this run (unlike the prior run, which surfaced `TypeError: Cannot read properties of null (reading 'length')` and `Input not found … "rules"`). |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 24s |
| grok-browser execution (scenario steps) | 12s |
| Execute via grok-browser (total) | 36s |
| Spec file generation | 25s |
| Spec script execution | 22s |
| **Total scenario run (with model)** | 2m 30s |

All rows are full-phase wall-clock (incl. model thinking and retries), not just tool latency.
The two `scenario steps` rows sum to `Execute via grok-browser (total)`. Reported Playwright runtime was 19.0s; the 22s figure includes `npx`/Node startup overhead.

## Summary

All three scenario steps passed in the MCP browser run against dev.datagrok.ai: `HELM.csv` opened (540 rows, HELM → Macromolecule), Bio > Manage > Monomer Libraries opened a "Manage Monomer Libraries" View (not a modal), 5 library checkboxes are shown (all checked by default), and toggling `HELMCoreLibrary123456.json` fired the "Monomer library user settings saved" balloon while shrinking the duplicate-symbols panel from 9 → 3 cards. The generated Playwright spec also PASSED (19.0s reported runtime). Login worked first try with the unmodified `loginToDatagrok` helper — no `<form>`-submit workaround was needed today, even though the prior run reported one. **Total scenario run (with model)**: ~2m 30s.

## Retrospective

### What worked well
- Bio > Manage > Monomer Libraries menu path is stable and reachable via `[name="div-Bio---Manage---Monomer-Libraries"]` after dispatching `mouseenter`+`mouseover` on the Manage submenu.
- The Manage Monomer Libraries View renders both the library list and the "Manage Duplicate Monomer Symbols" panel.
- Checkbox toggle is observable end-to-end: `checked` state flips, the duplicate-symbols panel resizes (9 ↔ 3 cards), and the balloon "Monomer library user settings saved" fires.
- Re-toggling cleanly restores state, so the scenario is idempotent across runs.
- `loginToDatagrok` (unmodified) worked on dev today — `keyboard.press('Enter')` did submit the form.

### What did not work
- Nothing failed this run.
- One stale balloon "21" appears alongside the expected ones — looks like a small badge/counter rather than a real notification; harmless but noisy in the balloon stack.

### Suggestions for the platform
- "Monomer lib updated" message: missing separator — `Monomer lib updated: PEPTIDE 336RNA 383` should be `PEPTIDE 336, RNA 383` (regression repro from prior run; still present).
- Login page (still applies as defensive hardening): wrap the Login/Password inputs in a real `<form>` and add `name="button-Login"` to the button. Today Enter worked, but the prior run's failure suggests a race or layout-dependent issue — making submit explicit and the button addressable would remove the variability.
- Investigate the stray `"21"` balloon — either suppress it or give it a meaningful label.

### Suggestions for the scenario
- Step 1: rename "sample_Helm.csv" to "HELM.csv" (the actual dataset path is `System:AppData/Bio/samples/HELM.csv`).
- Step 2: clarify that "Monomer Libraries" opens a **View**, not a modal dialog — the current "A dialog opens." line is misleading.
- Step 3: spell out the expected observable behaviour — e.g., "Uncheck one library and verify (a) the 'Monomer library user settings saved' balloon appears, (b) the Manage Duplicate Monomer Symbols panel updates, and (c) the change is reverted on re-check."
- Add an explicit cleanup step: re-check any library toggled during the test so user settings aren't mutated between runs.
