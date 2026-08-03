# Files & Sharing in Diff Studio — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open pk.ivp file preview | 30s | PASS | PASSED | `DiffStudio:previewIvp` function invoked with `FileInfo` from `grok.dapi.files.list('System:AppData/DiffStudio/library/', false, 'pk.ivp')`; `Func.prepare({file}) + call() + grok.shell.addView(view)` opens view 'PK' (TableView) with 10 input hosts: begin, end, step, dose, count, depot, centr, rateconstant, clearance, centralvolume. URL path becomes `/file/System.AppData/DiffStudio/library/pk.ivp?params:begin=0&end=12&step=0.01&dose=10000&count=1&…` |
| 2 | Set Step=0.1 (slider+textbox), Count=4 (clicker) | 30s | PASS | PASSED | Step via keyboard: click, Ctrl+A, "0.1", Tab → `0.1`. Count via `[name="input-host-count"] [name="icon-plus"]` clicked 3× (1→4). URL updates live to `step=0.1&count=4`. Playwright uses `{force: true}` on clicker because icons are hover-only |
| 3 | Open URL in new tab — same model & inputs | 30s | PASS | PASSED | `context.newPage()` + `goto(url)` loads same view 'PK' with `step=0.1` and `count=4` applied from URL params. Session cookie persists across tabs so login is not re-required |
| 4 | REMARK: with ≤3 curves, linechart only (no Multiaxis/Facet) | 20s | PASS | PASSED | Even at count=4, pk.ivp shows only 'PK' + 'Grid' tabs (no Multiaxis/Facet) — pk.ivp has 2 outputs (depot, centr) and the view apparently uses a count threshold different from Bioreactor/PK-PD. With count=1 confirmed: tabs are just PK + Grid, no Multiaxis or Facet — consistent with the REMARK |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 20s |
| grok-browser execution (scenario steps) | 1m 11s |
| Execute via grok-browser (total) | 2m 31s |
| Spec file generation | 1m 32s |
| Spec script execution | 1m 17s |
| **Total scenario run (with model)** | 5m 20s |

## Summary

Files & Sharing scenario reproduces fully on dev.datagrok.ai. pk.ivp loads via the `DiffStudio:previewIvp` function with a FileInfo from `grok.dapi.files.list('…/pk.ivp')`, Step edits via keyboard and Count clickers update the URL live, and reopening the URL in a fresh tab restores the same state (step=0.1, count=4). The REMARK is also verified: pk.ivp shows only PK + Grid tabs — no Multiaxis/Facet regardless of count. Playwright spec PASSED cleanly in 1.2m. **Total scenario run (with model)**: 5m 20s.

## Retrospective

### What worked well
- `DiffStudio:previewIvp` is a first-class function for opening .ivp files — scenario "click on the file pk.ivp" maps cleanly to `previewFn.prepare({file}).call()`
- URL round-trip via new tab validates the sharing contract: params after `?params:` are authoritative and apply on load
- Playwright `click({force: true})` bypasses the hover-visibility gate on `[name="icon-plus"]`/`[name="icon-minus"]` (learned from cyclic-models scenario — this is the right pattern for clicker-heavy flows)
- `context.newPage()` preserves the session cookie — no re-login needed

### What did not work
- pk.ivp does NOT render Multiaxis/Facet tabs at any count value tested (1, 4) — inconsistent with Bioreactor/PK-PD which switch to Multiaxis/Facet at count>=4. Either pk.ivp output shape (depot, centr) is below the threshold or the threshold is different
- Opening the file via DOM browse tree navigation (scenario wording "Browse > Files > App Data > Diff Studio > Library > pk.ivp") is not straightforward in automation — function-based preview path is used instead

### Suggestions for the platform
- Document the curve-count threshold for Multiaxis/Facet tabs per model; currently identical user action ("change count") produces different UI outcomes across Bioreactor vs pk.ivp
- Expose a `grok.shell.openFile(path)` API that mirrors the Browse tree click — this would let automation emulate the user action literally
- Session cookie for dev.datagrok.ai must be set on `*.datagrok.ai` (not just `dev.datagrok.ai`) so new tabs don't re-prompt for login — observed in Playwright where the new-tab login flow works only because same-context is used

### Suggestions for the scenario
- Step 1 should name the function that opens the preview (`DiffStudio:previewIvp`) or clarify what "click on the file" does if the file appears in the Browse tree
- Step 4 REMARK should specify for which model the threshold is ≤3 vs >=4 curves; currently the REMARK is ambiguous about whether it applies to pk.ivp or only Bioreactor/PK-PD
- Step 3 "Paste the copied URL into the address bar" is OS-specific; replace with "open a new browser tab and navigate to the copied URL"
