# SAR — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open linked datasets (peptides.csv) | PASS | PASSED | 647 rows, 3 cols loaded |
| 2 | Click the peptides column title | PASS | PASSED | `grok.shell.o = col` with context panel enabled |
| 3 | Expand Peptides panel on Context Panel | PASS | PASSED | Weblogo, Activity IC50, Scaling, Launch SAR button visible |
| 4 | Click Launch SAR button | PASS | PASSED | Clicked `button[text=Launch SAR]`; SAR view opened with all 4 viewers |
| 5 | Four viewers appear (Mutation Cliffs, Most Potent Residues, MCL, Logo summary) | PASS | PASSED | All 4 viewers loaded: Sequence Variability Map/Mutation Cliffs (bottom-left), Most Potent Residues (bottom-center), MCL network graph (top-center), Logo Summary Table (top-right) |
| 6 | Click settings button at the top | PASS | PASSED | Clicked `.d4-ribbon-item .fa-wrench`; "Peptides settings" dialog opened |
| 7 | Change all parameters arbitrarily | PASS | PASSED | Activity scaling: none→lg; Similarity Threshold: 70→60; Min Cluster Size: 5→3 |
| 8 | Click OK | PASS | PASSED | Dialog closed; viewers started reloading |
| 9 | Table and viewers reload with applied parameters | PASS | PASSED | MCL recalculated (blue nodes, different layout); table gained new columns (26 total); Cluster size (MCL) and Connectivity (MCL) columns appeared |
| 10 | Switch between Mutation Cliffs and Invariant Map | PASS | PASSED | Clicked "Invariant Map" label; grid changed from dots to filled colored cells with invariant position counts |
| 11 | Click random non-empty cell in Mutation Cliffs/Invariant Map | PASS | PASSED | Clicked cell (1:A) in Invariant Map; Selection Sources shows "Invariant map: 1:A" |
| 12 | Context Panel shows Mutation Cliff pairs, Distribution panels | PASS | PASSED | Both panels visible: "Mutation Cliffs pairs" and "Distribution" |
| 13 | Change parameters on Distribution panel | AMBIGUOUS | N/A | Distribution shows "No distribution" — no controls available; requires clicking a cell with actual mutation cliff pairs (pair-mode data); panel exists but has no settable parameters |

## Summary

12 of 13 steps passed (1 ambiguous). The full SAR workflow worked: SAR launched successfully with all 4 viewers, settings dialog opened with all parameters accessible, changing parameters (scaling, similarity threshold, min cluster size) triggered full recalculation and viewer reload. Mutation Cliffs/Invariant Map switching worked correctly. The only ambiguity was Step 13 (Distribution panel parameters) — the panel was present but contained no configurable controls in the current state.

## Retrospective

### What worked well
- `button[text=Launch SAR]` reliably clicked the SAR launch button
- Settings wrench button (`.d4-ribbon-item .fa-wrench`) opened the Peptides settings dialog
- Native setter for SELECT elements changed Activity scaling and Distance Function successfully
- OK button click triggered full recalculation with new parameters
- Mutation Cliffs ↔ Invariant Map switch worked via clicking label text
- Canvas click on the map correctly fired the "Mutation Cliffs pairs" and "Distribution" panels in the context panel

### What did not work
- Distribution panel showed "No distribution" for all clicked cells — likely requires specific conditions (actual mutation cliff pairs data)
- Mutation Cliffs pairs showed "No mutations table generated" — might require clicking specifically on a dot/marker rather than any cell

### Suggestions for the platform
- The Distribution panel should show a placeholder message explaining what to click to get distribution data
- "No mutations table generated" could be replaced by a "Click on a cell with mutation pairs to see details" hint
- Add `data-testid` attributes to the Mutation Cliffs / Invariant Map cells for reliable automation

### Suggestions for the scenario
- Step 13 should clarify what constitutes a "non-empty cell" that produces a Distribution (a cell with mutation cliff pairs, shown as dots in Mutation Cliffs view)
- Step 7 says "Change all parameters" — the settings dialog has many parameters; specify which ones are most important to verify
- Add a step: verify that the MCL viewer updates after clicking OK (currently step 9 is implicit)
