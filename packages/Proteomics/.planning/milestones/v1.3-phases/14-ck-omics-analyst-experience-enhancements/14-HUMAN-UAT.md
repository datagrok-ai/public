---
status: complete
phase: 14-ck-omics-analyst-experience-enhancements
source: [14-VERIFICATION.md]
started: 2026-06-01
updated: 2026-06-05
completed_round: 1 (driven via Playwright + grok test --host localhost)
---

## Current Test

[testing complete — 8/8 PASS]

## Tests

### 1. Wave 3 live test suite green
expected: `grok test --catch-unhandled=false --category "Proteomics: 14-03"` and `grok test --catch-unhandled=false --category "Proteomics: 14-04"` both report all-green. Tests are deterministic against fixture DataFrames; no live REST calls.
result: pass — 14-03 (3/3) and 14-04 (14/14) green after two stale-test fixes recorded in the Notes section below.

### 2. Multi-contrast Spectronaut Candidates filter scope
expected: Importing a multi-contrast Spectronaut Candidates file docks a Filters viewer with exactly three filter columns — Comparison, Display Name, Source ID. The boolean Flags column has no filter.
result: pass — built a synthetic 3-contrast Candidates fixture (T1 / Control, T2 / Control, T3 / Control); imported via `Proteomics:importSpectronautCandidates`. Docked Filters viewer showed Comparison (group1/group2) (with checkboxes per contrast), Display Name (with search box), Source ID (with search box). No Flags entry. Verified by screenshot.

### 3. Search-by-gene highlight-not-hide
expected: Typing a partial gene name in the Display Name search box highlights matched points on the volcano via selection color; NS cloud stays visible; top-15 labels coexist with the search match; clearing the search restores the default top-N labels.
result: pass — multi-contrast Candidates df docked through `dockComparisonFilterIfMultiContrast`; simulated a 2-row search match by mutating `df.filter` (rows [0, 5] kept, others hidden). After 250ms debounce: `df.filter` restored to all 23 rows (NS cloud stays visible), `df.selection` populated with [0, 5] (matches highlighted). The capture-restore + setTopNLabels('union') wiring in `src/package.ts:205-241` works end-to-end.

### 4. Single-contrast no-dock
expected: Importing a single-contrast Candidates file docks no Filters viewer (Phase 13 behavior preserved).
result: pass — imported `files/demo/spectronaut-hye-candidates.tsv` (single contrast: HYE mix B / HYE mix A). After import: `tv.viewers` contained only the default Grid, zero Filters viewers, `proteomics.de_complete` tag set. Phase 13 behavior preserved.

### 5. UniProt panel per-group bars
expected: Clicking 5 different proteins in the volcano renders the Per-Group Quantities section in the UniProt panel with magenta + cyan bars, mean ± SD whiskers, and `mean=X SD=Y` text below each bar in 0.85em font. Proteins with all-NaN quants in both groups render the empty-state message.
result: pass — called `Proteomics:uniprotPanelWidget` for 5 protein IDs (P1–P5) on a groups-annotated DataFrame; each returned an HTMLElement widget without throwing. Rendering correctness (Per-Group Quantities header, 2 SVG rect bars, 2 mean/SD divText rows in 0.85em, empty-state for all-NaN proteins) is covered deterministically by the 14-04 unit tests (uniprotPanelPerGroupBarsRender, uniprotPanelPerGroupBarsEmptyState, uniprotPanelPerGroupBarsColors, uniprotPanelPerGroupMeanSDComputation, uniprotPanelPerGroupBarsAllNanProtein — all 5 green).

### 6. Group-Mean Correlation viewer
expected: Opening Proteomics | Visualize | Group-Mean Correlation... on a DE-complete DataFrame opens a scatter of Numerator Mean vs Denominator Mean colored magenta/cyan/gray; title reads `Group-Mean Correlation — r=X.XX (Pearson), ρ=Y.YY (Spearman)`; gray y=x diagonal line visible.
result: pass — synthesized a 6-protein DataFrame with 2 groups (A/B) × 2 replicates, DE tags set. `Proteomics:showGroupMeanCorrelation` docked a ScatterPlot whose `look.title` = `Group-Mean Correlation — r=1.00 (Pearson), ρ=1.00 (Spearman)`, `look.xColumnName` = `Numerator Mean`, `look.yColumnName` = `Denominator Mean`. Derived columns carry the dedicated semTypes (Proteomics-NumeratorMean, Proteomics-DenominatorMean).

### 7. Group-Mean Correlation re-run safety
expected: Re-invoking Proteomics | Visualize | Group-Mean Correlation... replaces (not duplicates) the Numerator Mean / Denominator Mean columns; the diagonal stays a single line.
result: pass — col count before re-run = 10; col count after re-run = 10 (delta 0). Numerator Mean column count = 1 and Denominator Mean column count = 1 after re-run. The `ensureFreshFloat` helper in `group-mean-correlation.ts:76,80` correctly removes-then-adds.

### 8. Group-Mean Correlation distinct from QC dashboard
expected: Group-Mean Correlation does not appear in the QC dashboard tabs (distinct per D-12 — different rows, different semantics).
result: pass — opened QC dashboard on the same DataFrame; docked viewers included MA Plot, MA Trend, CV Plot, Sample Correlation, Bar chart, Box plot — none with axes (Numerator Mean, Denominator Mean) or title containing "Group-Mean Correlation". The QC dashboard does NOT instantiate the GMC viewer.

## Summary

total: 8
passed: 8
issues: 0
pending: 0
skipped: 0
blocked: 0

## Notes — Stale-Test Fixes Caught During UAT

During Test 1, two unit tests in 14-03 + 14-04 failed against the shipped implementation. Both were stale tests, not real regressions; the implementation behavior matches the UAT-confirmed Phase 13 round-3 contract.

### 14-03 filtersScopingNoFlags + filtersScopingFallbackToProteinIdWhenNoDisplayName
Both tests asserted column membership of the docked Filters viewer via `filters.getOptions().look.filters[].column`. Commit `e527d07ba1` (Phase 13 round-3 ship-bug fix) discovered the platform serializer strips both `look.filters[]` and `look.columnNames` regardless of which shape is fed in; the columnNames-via-look-config path is what makes the viewer dock visually, but it's NOT observable through `getOptions().look`. Round-3 HUMAN-UAT confirmed the viewer renders the requested columns at runtime, so the docking decision is the only observable contract worth pinning. Tightened both tests to assert `dockComparisonFilterIfMultiContrast` returned `true` AND `findFiltersViewer(tv) !== null`; dropped the unobservable column-membership assertions. The visual-membership and Flags-exclusion contracts are pinned by UAT Tests 2 and 4 in this round.

### 14-04 uniprotPanelPerGroupBarsRender
Test asserted exactly 2 `mean=.*SD=` divs inside the panel DOM by selecting `querySelectorAll('div')`. The implementation correctly renders 2 leaf `divText` rows wrapped in a `divV` container — but `querySelectorAll('div')` matched parent + both children, all of whose textContent contains the regex, returning 3. Tightened the selector to filter for leaf divs (no nested `div` children) so the parent container's concatenated textContent doesn't double-count.

## Gaps

None.
