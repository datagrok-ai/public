# Cloned Views — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Dataset**: System:AppData/Chem/tests/spgi-100.csv (100 rows, 88 columns)
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open spgi-100 | 6s | PASS | PASSED | `grok.dapi.files.readCsv` loaded 100 rows, 88 cols; 7 Molecule semtypes detected |
| 3 | Open the Filter Panel | 1s | PASS | PASSED | `tv.getFiltersGroup()`, 42 filter cards rendered |
| 4 | Navigate to Competition assay filter | <1s | PASS | PASSED | Histogram card present |
| 5 | Filter out missing values for Competition assay | 2s | PASS | PASSED | `fg.updateOrAdd({filterOutMissingValues: true})` → 80/100 rows |
| 6 | Set Stereo Category to S_ACHIR | 2s | PASS | PASSED | `fg.updateOrAdd({type: CATEGORICAL, selected:['S_ACHIR']})` → 28 rows |
| 7 | Sketch c1ccccc1 structure | 10s | PASS | PASSED | Clicked `.sketch-link` on Structure card, typed SMILES, Enter, OK → 11 rows |
| 8 | Clone View via View > Layout > Clone View | 7s | PASS | PASSED | `dispatchEvent(mouseenter)` on Layout submenu, clicked Clone View → "Table copy" |
| 9 | Verify cloned view filter state matches original | <1s | PASS | PASSED | 11 rows preserved, filter panel open, sketcher canvas present |
| 10 | Turn all filters off | 2s | PASS | PASSED | `.d4-filter-group-header input[type=checkbox]` click → 100 rows visible |
| 11 | Turn filters back on | 2s | PASS | PASSED | Re-click checkbox → 11 rows filtered again |
| 12 | Clear Structure, set to C1CCCCC1 | 10s | PASS | PASSED | `.chem-clear-sketcher-button` → 28; then C1CCCCC1 via sketch-link → 1 row |
| 13 | Remove Structure filter — state should not change | 2s | PASS | PASSED | `[name="icon-times"]` on Structure card, count stayed at 1 |
| 14 | Save the layout | 2s | PASS | PASSED | `grok.dapi.layouts.save()` → id `fbbf96e0-…` |
| 15 | Close the Filter Panel | 2s | PASS | PASSED | `fg.close()`, panel gone from DOM |
| 16 | Apply saved layout — panel opens without Structure | 5s | PASS | PASSED | `loadLayout(saved)` → filter panel open, no Structure, 1 row preserved |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m 28s |
| grok-browser execution (scenario steps) | 48s |
| Execute via grok-browser (total) | 3m 16s |
| Spec file generation | 14s |
| Spec script execution | 55s |
| **Total scenario run (with model)** | 5m 44s |

Note: `Spec file generation` is review-only — the existing `cloned-views-spec.ts` was validated against the run log but not rewritten per the invoker's request. `Spec script execution` is the wall-clock of `npx playwright test cloned-views-spec.ts --config playwright.config.files-and-sharing.ts --headed` (test itself: 51.8s; Playwright startup + report: ~3s).

## Summary

All 15 scenario steps PASSed on dev. spgi-100.csv loads correctly this time (previous run had to substitute SPGI.csv). Cloned view correctly inherits every filter from the original (Competition assay missing-values, Stereo Category=S_ACHIR, Structure=c1ccccc1). Toggling filters off/on, re-sketching the structure, removing the Structure filter card, and layout save/close/apply all behave as specified. The existing Playwright spec replayed the scenario end-to-end in 51.8s without modification. Total scenario run (with model): **5m 44s**.

## Retrospective

### What worked well
- `grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv')` succeeded on dev — the file is now available (it wasn't in the previous 2026-04-06 run)
- Scoping the Structure sketch-link by filter header text avoided ambiguity across 7 Molecule columns
- `fg.updateOrAdd({filterOutMissingValues: true})` + `fg.updateOrAdd({selected: ['S_ACHIR']})` hit the filter state deterministically without UI submenu navigation
- Clone View via `dispatchEvent(mouseenter)` on the Layout submenu node remains the reliable way to trigger the nested menu path
- Existing `cloned-views-spec.ts` replayed cleanly against dev — no regressions

### What did not work
- Running `npx playwright test` without `--config playwright.config.files-and-sharing.ts` fails with "No tests found" because Playwright's default testMatch requires `.spec.ts` / `.test.ts` (with a dot) but these files use `-spec.ts` (with a hyphen). The spec is not discoverable without the root config
- Chrome DevTools MCP `take_snapshot` fires a reference-file reminder hook even when the references were already read in this session — not a blocker but a friction point

### Suggestions for the platform
- Add `name=` attribute to the filter panel title-bar close icon so it's reachable via selector (currently only `fg.close()` works reliably)
- Consider a JS API like `fg.updateOrAdd({type: 'Chem:substructureFilter', column: 'Structure', smiles: 'c1ccccc1'})` to avoid having to drive the sketcher dialog for substructure filters in automation

### Suggestions for the scenario
- Step numbering jumps from 1 to 3 (step 2 is missing) — renumber or mark step 2 explicitly as a no-op
- Clarify step 13: "filtered state should not change" — the count stays the same because the original view still holds the Structure filter on the shared DataFrame. Spelling this out prevents tester confusion
- Pre-condition that spgi-100.csv must exist on the target server — if it's missing, scenario cannot run faithfully (2026-04-06 run had to substitute SPGI.csv; 2026-04-22 run found the file available)
