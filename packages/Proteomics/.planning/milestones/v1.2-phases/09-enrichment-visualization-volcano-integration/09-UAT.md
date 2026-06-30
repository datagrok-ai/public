---
status: complete
phase: 09-enrichment-visualization-volcano-integration
source: [09-01-SUMMARY.md, 09-02-SUMMARY.md]
started: 2026-03-07T02:00:00Z
updated: 2026-03-07T02:05:00Z
---

## Current Test

[testing complete]

## Tests

### 1. Enrichment Auto-Opens After Analysis
expected: After running enrichment analysis (Proteomics > Analyze > Enrichment Analysis), clicking OK completes the g:Profiler query and automatically opens enrichment dot plot and bar chart viewers alongside the enrichment results table. No manual step needed to see visualizations.
result: pass

### 2. Enrichment Dot Plot
expected: The dot plot shows enriched terms on the Y axis and a numeric measure (e.g., -log10 FDR) on the X axis. Dot size reflects gene count (or intersection size) and dot color reflects significance. Only the top N terms are shown (not all hundreds).
result: pass

### 3. Enrichment Bar Chart
expected: The bar chart shows top enriched terms as horizontal bars sorted by significance (most significant at top). Term names are readable and truncated if too long.
result: pass

### 4. Cross-DataFrame Volcano Wiring
expected: When you click on an enrichment term row in the enrichment results table, the corresponding member genes/proteins are highlighted (selected) in the volcano plot of the original protein DataFrame. This links enrichment results back to the DE data.
result: pass

### 5. Re-Open Menu Entry
expected: With an enrichment results table active, the menu entry "Proteomics > Visualize > Enrichment Charts..." is available. Clicking it re-opens the dot plot and bar chart viewers for that enrichment table, with volcano wiring if the protein DE table is still open.
result: pass

## Summary

total: 5
passed: 5
issues: 0
pending: 0
skipped: 0

## Gaps

[none]
