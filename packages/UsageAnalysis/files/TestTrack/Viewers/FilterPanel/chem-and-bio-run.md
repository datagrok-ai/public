# Chem and Bio — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

### Section 1: Chem filter

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open spgi-100 | 10s | PASS | PASSED | `grok.dapi.files.readCsv` + `addTableView`; semType `Molecule` detected on 7 columns; canvas + 5s settle |
| 2 | Open the Filter Panel | 22s | PASS | PASSED | UI click on `[name="div-section--Filters"]` did not open panel — JS API fallback `tv.getFiltersGroup()` opened 42 cards including Structure |
| 3 | Draw CCC(N(C)C)=O in Structure filter | 7s | PASS | PASSED | `.sketch-link` click → SMILES typed via native value setter + Enter → `[name="button-OK"]`; 14 rows |
| 4 | Click gear icon, reveal search type dropdown | 5s | PASS | PASSED | `.chem-search-options-icon` inside Structure card; `<select>` with 7 options revealed |
| 5a | Switch to Included in | 4s | PASS | PASSED | 0 rows asserted |
| 5b | Switch to Exact | 4s | PASS | PASSED | 0 rows asserted |
| 5c | Switch to Similar | 4s | PASS | PASSED | 0 rows asserted |
| 5d | Switch to Not contains | 4s | PASS | PASSED | 86 rows asserted |
| 5e | Switch to Not included in | 4s | PASS | PASSED | 100 rows asserted |
| 5f | Switch to Contains | 4s | PASS | PASSED | 14 rows asserted |
| 6 | Close All | 3s | PASS | PASSED | `grok.shell.closeAll()` |

### Section 2: Bio

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open peptides.csv | 7s | PASS | PASSED | 647 rows × 3 cols; `AlignedSequence:Macromolecule` detected |
| 2 | Wait for AlignedSequence to render | (incl. step 1) | PASS | PASSED | Canvas poll + 5s settle baked into open phase |
| 3 | Open the Filter Panel | 16s | PASS | PASSED | `tv.getFiltersGroup()`; AlignedSequence card present on first try (no close/reopen needed this run) |
| 4 | Enter T-T-Y-K-N-Y-V in Substructure | 18s | PASS | PASSED | `input[placeholder="Substructure"]`, value via native setter + Enter; 28 rows |
| 5 | Close All | 7s | PASS | PASSED | `grok.shell.closeAll()` |

**Time** = step 2b wall-clock per step (incl. model thinking + waits). **Result** = step 2b outcome. **Playwright** = step 2e outcome (existing spec ran without modification).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 52s |
| grok-browser execution (scenario steps) | 58s |
| Execute via grok-browser (total) | 2m 50s |
| Spec file generation | 42s |
| Spec script execution | 47s |
| **Total scenario run (with model)** | 5m 16s |

`Execute via grok-browser (total)` = end of step 2b – start of step 2b. The two `scenario steps` rows (thinking + grok-browser execution) sum to that total. `Spec script execution` = `npx playwright test ... --headed` wall-clock (44.9s test + ~3s startup/teardown). `Total scenario run` = end of step 2e – start of step 2b.

## Summary

Ran chem-and-bio scenario end-to-end against dev. All 11 scenario steps passed in the MCP browser phase (Chem: open spgi-100 → filter panel → sketch CCC(N(C)C)=O → 14 rows → gear → 6 search types matched 0/0/0/86/100/14 → close. Bio: open peptides → filter panel → T-T-Y-K-N-Y-V → 28 rows → close). Existing `chem-and-bio-spec.ts` was left untouched per user instruction; re-running it produced 7 PASSED softSteps in 44.9s. **Total scenario run (with model)**: 5m 16s.

## Retrospective

### What worked well
- All chem search-type counts matched the scenario exactly (0/0/0/86/100/14)
- Bio filter panel opened on the first `getFiltersGroup()` call this run — close/reopen fallback in the spec was not triggered
- Native value setter (`Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set`) reliably triggered the Dart change listener for both SMILES and Substructure inputs
- Existing spec was already correct against dev — no regeneration needed; running it gives a sub-minute regression check

### What did not work
- UI click on `[name="div-section--Filters"]` returned `'clicked'` but did not open `[name="viewer-Filters"]` — required `grok.shell.tv.getFiltersGroup()` JS API fallback. Spec uses this fallback directly, sidestepping the issue
- `npx playwright test <path>` requires `--config=playwright.config.files-and-sharing.ts` because the repo's specs use the `*-spec.ts` naming pattern (hyphen) which doesn't match Playwright's default `*.spec.ts` (dot) pattern

### Suggestions for the platform
- `[name="div-section--Filters"]` click should open the Filters panel reliably (or the toolbox click handler should be inspected — it currently appears to no-op for the Filters section in some sequences)

### Suggestions for the scenario
- Step 4 should clarify that "click the gear icon" reveals an inline `<select>` dropdown directly on the filter card (not a separate dropdown popup); the current wording is accurate but mildly ambiguous
- Step 5 ordering (Included in / Exact / Similar / Not contains / Not included in / Contains) does not match the dropdown's native option order (Contains / Included in / Exact / Stereo agnostic / Similar / Not contains / Not included in) — the test still works, but listing values in dropdown order would be slightly more natural for manual execution
