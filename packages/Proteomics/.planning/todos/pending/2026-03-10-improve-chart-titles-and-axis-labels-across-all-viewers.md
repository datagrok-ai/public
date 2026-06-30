---
created: 2026-03-10T14:16:05.170Z
title: Improve chart titles and axis labels across all viewers
area: ui
files:
  - packages/Proteomics/src/viewers/volcano.ts
  - packages/Proteomics/src/viewers/pca-plot.ts
  - packages/Proteomics/src/viewers/heatmap.ts
  - packages/Proteomics/src/viewers/enrichment-viewers.ts
  - packages/Proteomics/src/viewers/qc-dashboard.ts
---

## Problem

Charts across the Proteomics package lack descriptive axis labels and could have better titles. While Phase 11 added viewer titles to volcano, PCA, and heatmap, the axis labels on all viewers (including QC dashboard charts and enrichment visualizations) still use raw column names or default Datagrok labels rather than human-readable descriptions (e.g., "log2(Fold Change)" instead of "log2FC", "Sample" instead of column key).

## Solution

Audit all viewer creation calls across the package and set meaningful:
- **X/Y axis labels** with proper units and formatting (e.g., "-log10(p-value)", "log2(Fold Change)", "PC1 (X% variance)")
- **Titles** for any viewers not yet titled (enrichment dot plots, QC dashboard charts)
- Use Datagrok viewer `setOptions()` for axis label properties where available
