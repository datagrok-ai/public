# R-Groups Analysis — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open sar_small.csv (200 rows, `smiles` semType=Molecule) | 8s | PASS | PASSED | Grid renders with canvas |
| 2 | Chem → Analyze → R-Groups Analysis → dialog opens | 2s | PASS | PASSED | Dialog title "R-Groups Analysis" |
| 3 | Click MCS → scaffold computed into sketcher | 8s | PASS | PASSED | No error balloon |
| 4 | Click OK → Trellis plot + R1/R2/R3/R4 columns | 12s | PASS | PASSED | Viewers: Grid, Trellis plot; new cols: Core, R1-R4, r-groups-highlight_0, isMatch |
| 5 | Open smiles.csv → R-Groups → MCS → OK → no trellis | 10s | PASS | PASSED | Viewers: only Grid (no R Groups found, no trellis produced) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 30s |
| grok-browser execution (scenario steps) | 50s |
| Execute via grok-browser (total) | 2m 20s |
| Spec file generation | 45s |
| Spec script execution | 58.7s |
| **Total scenario run (with model)** | ~4m |

## Summary

R-Groups Analysis works on sar_small.csv: MCS auto-populates the sketcher, OK produces a Trellis plot and appends R1–R4 / Core columns. On smiles.csv (no common scaffold), the operation correctly produces no Trellis plot. Playwright replay passes both cases in under a minute.

## Retrospective

### What worked well
- Chem menu items reachable by name attr (`div-Chem`) + submenu text match on `.d4-menu-item-label`
- MCS button found by text (`button:has-text("MCS")`) inside `.d4-dialog`
- 8s wait between MCS click and OK is sufficient for 200-row sar_small
- Trellis plot / R-group column presence are both strong assertions

### What did not work
- `grok.shell.warnings` / `grok.shell.infos` stay empty even when the "No R Groups were found" balloon fires — balloons self-expire before the JS scan catches them. Alternative: query `.grok-balloon` DOM immediately after OK click

### Suggestions for the platform
- Persist recent balloons on `grok.shell.warnings` / `grok.shell.balloons` so automation can assert "expected balloon appeared" after the visual fade-out
- Expose an explicit JS API for R-Groups (`Chem:rGroupsAnalysis(molCol, core)`) so automation doesn't have to drive the dialog

### Suggestions for the scenario
- Step 11 ("without clicking MCS") → add explicit note that the "No core was provided" balloon also fades quickly; use grok.shell warnings log check instead
- Clarify: Replace latest = **checked** ⇒ new results REPLACE existing columns; unchecked ⇒ results are APPENDED (the wording currently leads with the negative case)
- Specify "trellis plot" as the expected viewer type so QA can assert it verbatim
