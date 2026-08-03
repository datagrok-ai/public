<!-- cSpell:words FASTA UMAP Datagrok datagrok dataframe Macromolecule Hamming Levenshtein Needlemann Wunsch dapi DBSCAN -->
# Sequence Space — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open sample_FASTA.csv | 6s | PASS | PASSED | Loaded as `System:AppData/Bio/samples/FASTA.csv` (actual filename `FASTA.csv`, not `sample_FASTA.csv`); 64 rows; `Sequence` detected as Macromolecule. HELM.csv (540 rows, `HELM` col) and MSA.csv (540 rows, `MSA` col) also loaded cleanly |
| 2 | Open Bio > Search > Sequence Space | 5s | PASS | PASSED | Dialog "Sequence Space" opened. Menu path on the live platform is `Bio > Analyze > Sequence Space...` (not `Bio > Search`). Used the actual path; intent achieved |
| 3 | Click OK to run with default parameters | 5s | PASS | PASSED | Defaults: Column auto-selected (Sequence/HELM/MSA), Preprocessing=Encode Sequences, Method=UMAP, Similarity=Hamming, Plot embeddings=true, Cluster embeddings=true. Scatter plot added; +3 columns (`Embed_X_1`, `Embed_Y_1`, `Cluster (DBSCAN)`) |
| 4 | Re-open Bio > Search > Sequence Space | 3s | PASS | PASSED | Dialog re-opened. Same menu-path remap as step 2 |
| 5 | Change Similarity and Method arbitrarily | <1s | PASS | PASSED | Method UMAP → t-SNE, Similarity Hamming → Levenshtein (Needlemann-Wunsch on FASTA only) via `<select>` + `change` event |
| 6 | Click OK to run with edited parameters | 3s | PASS | PASSED | Second scatter plot rendered; +3 more columns (`Embed_X_2`, `Embed_Y_2`, `Cluster (DBSCAN) (2)`); no console warnings |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~2m |
| grok-browser execution (scenario steps) | ~50s |
| Execute via grok-browser (total) | ~2m 50s |
| Spec file generation | ~30s |
| Spec script execution | 1m 0s (FASTA 19.9s + HELM 22.1s + MSA 19.0s) |
| **Total scenario run (with model)** | ~4m 20s |

## Summary

All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets: FASTA, HELM, MSA). Sequence Space reliably produces a Scatter plot of the embedding space and adds `Embed_X_*`, `Embed_Y_*`, and `Cluster (DBSCAN)*` columns to the dataframe on every run. The scenario's menu-path ("Search" → actual "Analyze") and file-name (`sample_*.csv` → actual `*.csv`) discrepancies are wording issues, not test failures. **Total scenario run (with model)** ≈ 4m 20s.

## Retrospective

### What worked well
- JS API dataset open via `grok.dapi.files.readCsv('System:AppData/Bio/samples/<FASTA|HELM|MSA>.csv')` worked first try; Bio package initialised in time after the 5s settle.
- The Sequence Space dialog auto-selected the Macromolecule column (`Sequence` / `HELM` / `MSA`) for each dataset — no user input required.
- `<select>` dropdowns inside the dialog accepted programmatic value changes when the `change` event was dispatched — no keyboard simulation required.
- Playwright `page.locator(...).hover()` correctly opens top-menu submenus, so the spec did not need to fall back to JS API for menu navigation (unlike the earlier MCP run where only `chrome-devtools__hover` triggered the Dart `onMouseEnter` listener).
- All three notations (FASTA / HELM / MSA) produced an "Embeddings" scatter plot under both UMAP/Hamming defaults and t-SNE/Levenshtein edits, confirming notation-agnostic behaviour.
- Each run is strictly additive on the dataframe: default run adds `Embed_X_1`/`Embed_Y_1`/`Cluster (DBSCAN)`, edited run adds `Embed_X_2`/`Embed_Y_2`/`Cluster (DBSCAN) (2)` — handy for before/after comparisons.

### What did not work
- The scenario menu path "Bio > Search > Sequence Space" does not exist on the live platform — the function lives at "Bio > Analyze > Sequence Space..." (same discrepancy as the sibling Activity Cliffs scenario). I used the actual path; the step still PASSED because the intent was achieved.
- The scenario file names "sample_FASTA.csv / sample_HELM.csv / sample_MSA.csv" do not match the real files, which are simply `FASTA.csv` / `HELM.csv` / `MSA.csv` under `System:AppData/Bio/samples/`.
- CDP synthetic clicks on top-menu items did not open submenus during the MCP run; only `chrome-devtools__hover` (real mouse move) triggered the Dart `onMouseEnter` listener that calls `_showSubMenu`. Workaround for nested submenus: hover the parent via MCP, then dispatch `mouseenter`/`mouseover` to the child in `evaluate_script`.
- The dialog's auto-picked column is the first Macromolecule column, not necessarily the one the user wants if a dataframe has more than one — worth noting for scenarios with multi-sequence tables.

### Suggestions for the platform
- Consider making the d4 top menu open on click *or* hover. Right now CDP / synthetic clicks without hover do not open the submenu, which makes test automation outside Playwright harder than it needs to be.
- Add a `name=` attribute on the resulting "Embeddings" Scatter plot viewer (or its header) so spec assertions can target it specifically rather than just asserting "any Scatter plot exists".
- Sequence Space always adds `(DBSCAN)` clustering columns by default — that makes sense for exploration, but a `Cluster embeddings` toggle that is off by default (with the checkbox already present) would reduce dataframe pollution on repeated runs.
- The default Similarity is `Hamming`, which only really works when all sequences have equal length. On variable-length FASTA inputs this can silently degrade quality — consider picking `Levenshtein` or `Needlemann-Wunsch` by default when lengths differ.

### Suggestions for the scenario
- Update the menu path from "Bio > Search > Sequence Space" to the actual "Bio > Analyze > Sequence Space..." — the current text sends testers to a wrong submenu.
- Update the file names from "sample_FASTA.csv / sample_HELM.csv / sample_MSA.csv" to "FASTA.csv / HELM.csv / MSA.csv" (paths under `System:AppData/Bio/samples/`).
- Specify which Similarity / Method values to switch to (e.g., "Levenshtein" and "t-SNE") rather than "arbitrarily" — that would make the scenario reproducible across runs and give a clean before/after for the embedding comparison.
- Add a verification step after each OK click ("observe that a Scatter plot appears and 3 new columns — Embed_X, Embed_Y, Cluster (DBSCAN) — are added to the dataframe"), so the tester knows explicitly what to check.
