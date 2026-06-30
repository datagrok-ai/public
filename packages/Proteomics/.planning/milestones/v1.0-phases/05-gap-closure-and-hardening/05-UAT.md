---
status: diagnosed
phase: 05-gap-closure-and-hardening
source: [05-01-SUMMARY.md, 05-02-SUMMARY.md, 05-03-SUMMARY.md]
started: 2026-03-03T19:45:00Z
updated: 2026-03-05T10:05:00Z
---

## Current Test

[testing complete]

## Tests

### 1. Semantic Type Auto-Detection on File Open
expected: Open a MaxQuant proteomics file (e.g. cptac-spike-in.txt). Intensity columns — including those with "Log2 LFQ intensity" prefix (space format) — should be auto-detected as Proteomics-Intensity semType. Protein ID columns should be detected as Proteomics-ProteinId. Check column properties to confirm.
result: pass

### 2. DE Analysis Executes Without Errors
expected: Run Proteomics > Analyze > Differential Expression on a loaded proteomics dataset. The analysis should complete and add fold-change and p-value columns to the table. No "function not found" errors from R script call name mismatches.
result: pass

### 3. Heatmap Creation Shows Progress Indicator
expected: After DE completes, run Proteomics > Visualize > Expression Heatmap. A "Creating heatmap..." progress indicator should appear in the Datagrok taskbar while the heatmap loads. The UI should not appear frozen.
result: pass

### 4. Heatmap Filter Isolation from Volcano Plot
expected: With both a heatmap and volcano plot open, the heatmap displays only top-N significant proteins while the volcano plot retains ALL data points. Filtering/selection in the heatmap should NOT reduce or affect the volcano plot's rows.
result: pass

### 5. Dendrogram Tree on Heatmap
expected: The expression heatmap should display a dendrogram (tree diagram) on the left side showing hierarchical clustering of proteins/samples. Rows should be reordered by cluster similarity. The tree should show branching lines (~150px wide).
result: issue
reported: "fail. no dendrogram appears. additionally for the cptac-spike demo data set no gene name is shown in the heatmap"
severity: major

### 6. DE Fallback When R/Limma Unavailable
expected: If limma R package is not available, the DE analysis should fall back to client-side t-test quickly (within seconds, not 30+ seconds). A warning should indicate the fallback occurred.
result: pass

## Summary

total: 6
passed: 5
issues: 1
pending: 0
skipped: 0

## Gaps

- truth: "Heatmap row labels show gene names (or protein IDs as fallback) for all datasets"
  status: failed
  reason: "User reported: for the cptac-spike demo data set no gene name is shown in the heatmap"
  severity: major
  test: 5
  root_cause: "cptac-spike-in.txt has no gene name column. heatmap.ts:118 findColumn returns null for GENE_SYMBOL and silently shows no row labels — no fallback to protein ID column."
  artifacts:
    - path: "packages/Proteomics/src/viewers/heatmap.ts"
      issue: "Lines 118-122: no fallback when gene symbol column not found"
    - path: "packages/Proteomics/files/demo/cptac-spike-in.txt"
      issue: "Demo file has no gene name column"
  missing:
    - "Add fallback to protein ID column when gene symbol not found in heatmap label logic"

- truth: "Heatmap displays dendrogram tree showing hierarchical clustering with branching lines"
  status: failed
  reason: "User reported: no dendrogram appears"
  severity: major
  test: 5
  root_cause: "Dendrogram rendering code (injectTreeForGrid) is entirely commented out at lines 150-168. The replacement dFunc.apply approach passes wrong argument types (Grid instead of DataFrame, string[] instead of ColumnList) and fails silently via catch block."
  artifacts:
    - path: "packages/Proteomics/src/viewers/heatmap.ts"
      issue: "Lines 140-143: wrong args (grid vs heatmapDf, string[] vs ColumnList). Lines 150-168: correct tree rendering code commented out."
  missing:
    - "Restore commented-out TreeHelper/DendrogramService approach or fix dFunc.apply argument types"
