<!-- cSpell:words FASTA UMAP SALI Datagrok datagrok dataframe Macromolecule Hamming Levenshtein Needlemann Wunsch dapi initialised normalised -->
# Sequence Activity Cliffs — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open sample_FASTA.csv | 6s | PASS | PASSED | Loaded as `System:AppData/Bio/samples/FASTA.csv` (actual filename `FASTA.csv`, not `sample_FASTA.csv`); 64 rows; Sequence detected as Macromolecule |
| 2 | Open Bio > Search > Sequence Activity Cliffs | 5s | PASS | PASSED | Dialog "Activity Cliffs" opened. Function is internally named "Sequence Activity Cliffs" but its menu path on the live platform is `Bio > Analyze > Activity Cliffs...` (not `Bio > Search`). Used the actual path; intent achieved |
| 3 | Click OK to run with default parameters | 6s | PASS | PASSED | Defaults: Column=Sequence, Activities=Length, Method=UMAP, Similarity=Hamming, cutoff=80. Scatter plot ("Activity cliffs", 2 cliffs) added; column count 6 → 9 (embedding X/Y + SALI) |
| 4 | Re-open Bio > Search > Sequence Activity Cliffs | 4s | PASS | PASSED | Dialog re-opened. Same menu-path remap as step 2 |
| 5 | Change Similarity and Method arbitrarily | <1s | PASS | PASSED | Method UMAP → t-SNE, Similarity Hamming → Levenshtein via `<select>` + `change` event |
| 6 | Click OK to run with edited parameters | 5s | PASS | PASSED | Second scatter plot rendered; total scatter viewer count = 2; column count 9 → 12; no console warnings |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~2m |
| grok-browser execution (scenario steps) | ~26s |
| Execute via grok-browser (total) | ~2m 26s |
| Spec file generation | ~1m |
| Spec script execution | 1m 0s (FASTA 18.3s + HELM 22.1s + MSA 19.8s) |
| **Total scenario run (with model)** | ~4m 26s |

## Summary

All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets). Activity Cliffs reliably produces an "Activity cliffs" Scatter plot and adds embedding + SALI columns to the dataframe. The two menu/file-name discrepancies in the scenario text (path "Search" → actual "Analyze"; filename `sample_FASTA.csv` → actual `FASTA.csv`) are scenario wording issues, not test failures — captured under "Suggestions for the scenario" below. **Total scenario run (with model)** ≈ 4m 26s.

## Retrospective

### What worked well
- JS API dataset open via `grok.dapi.files.readCsv('System:AppData/Bio/samples/FASTA.csv')` worked first try; Bio package initialised in time after the 5s settle.
- The Activity Cliffs dialog defaulted Column to the Macromolecule column (Sequence) automatically.
- `<select>` dropdowns inside the dialog accepted programmatic value changes when the `change` event was dispatched — no keyboard simulation required.
- Playwright `page.locator(...).hover()` correctly opens top-menu submenus, so the spec did not need to fall back to JS API for menu navigation.
- All three sample datasets (FASTA / HELM / MSA) produced an Activity-cliffs scatter plot, confirming notation-agnostic behavior.

### What did not work
- The scenario menu path "Bio > Search > Sequence Activity Cliffs" does not exist on the live platform — the function lives at "Bio > Analyze > Activity Cliffs..." (the function name is "Sequence Activity Cliffs", but its menu location is under Analyze). I used the actual path; the step still PASSED because the intent was achieved.
- The scenario file names "sample_FASTA.csv / sample_HELM.csv / sample_MSA.csv" do not match the real files, which are simply `FASTA.csv` / `HELM.csv` / `MSA.csv` under `System:AppData/Bio/samples/`.
- CDP synthetic clicks on top-menu items did not open submenus during the MCP run; only `chrome-devtools__hover` (which moves the real mouse) triggered the Dart `onMouseEnter` listener that calls `_showSubMenu`. Quirk of the d4 horizontal menu's hover-to-open behavior.
- The "Activities" input defaulted to the `Length` column (an int) rather than the `Activity` column, which is semantically wrong for an "Activity Cliffs" analysis. The default still ran without error, but the result is not biologically meaningful.

### Suggestions for the platform
- Improve default-column detection in the Activity Cliffs dialog: when a column literally named "Activity" exists and is numeric, prefer it over other numeric columns (currently picks the first numeric column, which here is `Length`).
- Consider making the d4 top menu open on click *or* hover. Right now CDP / synthetic click without hover never opens the submenu, which makes test automation outside Playwright harder than it needs to be.
- Add a `name=` attribute on the resulting "Activity cliffs" Scatter plot viewer header so spec assertions can target it more specifically than just "first Scatter plot".

### Suggestions for the scenario
- Update the menu path from "Bio > Search > Sequence Activity Cliffs" to the actual "Bio > Analyze > Activity Cliffs..." (and clarify that the function is internally registered as "Sequence Activity Cliffs").
- Update the file names from "sample_FASTA.csv / sample_HELM.csv / sample_MSA.csv" to "FASTA.csv / HELM.csv / MSA.csv" (paths under `System:AppData/Bio/samples/`).
- Add a step instructing the tester to set the Activities column to `Activity` before clicking OK on the default run, so the result reflects the intent of the test rather than the platform's default-column choice.
- Specify which Similarity / Method values to switch to (e.g., "Levenshtein" and "t-SNE") rather than "arbitrarily" — that would make the scenario reproducible across runs.
