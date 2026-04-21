# Sensitivity Analysis in Diff Studio (Bioreactor) — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open Diff Studio + Bioreactor; Edit toggle OFF | 30s | PASS | PASSED | `runDiffStudio` + `grok.shell.addView`, double-click `.diff-studio-hub-card` (Bioreactor), verified Edit toggle `ui-input-switch` NOT `ui-input-switch-on`; Process mode input host present |
| 2 | Click Sensitivity icon — SA view opens | 15s | PASS | PASSED | Dispatched mousedown/mouseup/click on `.diff-studio-ribbon-sa-icon`. View switches to 'Bioreactor - comparison' (TableView) with body text containing "Sensitivity Analysis", Method=Monte Carlo, Samples=10, Process mode and 21 input hosts visible, 43 `.ui-input-switch` switchers for parameter selection |
| 3 | Modify Process mode; FFox/KKox cascade | 15s | PASS | PASSED | `select.value='Mode 1'` + input/change events → FFox 0.2→0.163, KKox 0.2→0.24, MEAthiol 15→10, Gas 1→0.5 (4 inputs cascaded); parameters remain selectable throughout |
| 4 | Run SA; 4 viewers visualize results | 20s | PASS | PASSED | Enabled FFox and KKox switchers (`.ui-input-switch` click turns them on), clicked `.fa-play` on ribbon. After ~15s: 5 `.d4-viewer` elements — Correlation plot, Grid, PC plot, Scatter plot (scenario expects 4). Status bar shows "Columns: 34, Rows: 10" |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 23s |
| grok-browser execution (scenario steps) | 40s |
| Execute via grok-browser (total) | 2m 3s |
| Spec file generation | 44s |
| Spec script execution | 39s |
| **Total scenario run (with model)** | 3m 26s |

## Summary

Sensitivity Analysis scenario fully reproduces against dev.datagrok.ai. Bioreactor loads from the DiffStudio hub (library tile), clicking the `.diff-studio-ribbon-sa-icon` switches to the SA comparison view with Monte Carlo as default method, Process mode cascade updates 4 dependent inputs (FFox, KKox, MEAthiol, Gas), and Run produces the expected 4 visualizations (Correlation plot, Grid, PC plot, Scatter plot; 5 `.d4-viewer` total including a secondary Grid). Playwright spec PASSED cleanly in 35.2s. **Total scenario run (with model)**: 3m 26s.

## Retrospective

### What worked well
- `.diff-studio-ribbon-sa-icon` is a specific, stable CSS class on the Sensitivity icon — no text-matching needed
- Process mode cascade is consistent between Diff Studio (Bioreactor) and Sensitivity Analysis views — same `Mode 1` mapping produces same 4-input cascade
- `input.ui-input-editor` scoping on `[name="input-host-*"]` avoids the strict-mode violation with range sliders (learned from the scripting scenario)
- Run produces visible viewers within ~15s; a single `waitForTimeout(15000)` is sufficient at 10 samples

### What did not work
- Scenario step 3 wording "use the switchers to select the following parameters" conflates two different things: the Process mode dropdown and the per-input switchers on the left of each host. The cascade verified is via Process mode change only; the param switchers are a separate toggle for "include this input in the analysis"
- Run requires at least one input switcher turned ON (otherwise SA produces an empty grid); this is not stated in the scenario
- The scenario says "four viewers should open" but dev actually opens 5 `.d4-viewer` elements (Correlation, Grid, PC plot, Scatter plot — the Grid shows sample input/output values; a secondary Grid viewer contains the full SA results table)

### Suggestions for the platform
- Add `name=` attributes to SA ribbon buttons (Run, Stop, Method selector) so automation does not rely on Font-Awesome class discrimination
- Pre-enable a sensible default input switcher (e.g., FFox) when entering SA view so Run is immediately actionable without requiring the user to toggle inputs first
- Distinguish the "primary" 4 result viewers from the auxiliary Grid so the scenario count matches what users see

### Suggestions for the scenario
- Step 3 should separately describe (a) Process mode change cascades, and (b) the per-input switchers for SA inclusion — they are confusingly combined in the current wording
- Step 4 should state the Run precondition: "ensure at least one input switcher is ON before clicking Run"
- Clarify the viewer count — expected 4 primary visualizations (Correlation plot, PC plot, Scatter plot, Grid) plus an auxiliary result Grid (5 total)
