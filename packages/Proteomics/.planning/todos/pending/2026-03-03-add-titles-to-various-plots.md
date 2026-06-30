---
created: 2026-03-03T19:12:00.000Z
title: Add titles to various plots
area: ui
files:
  - packages/Proteomics/src/viewers/volcano.ts
  - packages/Proteomics/src/viewers/pca-plot.ts
  - packages/Proteomics/src/viewers/heatmap.ts
---

## Problem

The volcano plot, PCA scatter, and expression heatmap viewers are missing descriptive titles. Users need to quickly identify what each plot shows, especially when multiple viewers are open simultaneously.

## Solution

Add titles to each visualization:
- **Volcano plot**: e.g. "Volcano Plot: Treatment vs Control"
- **PCA plot**: e.g. "PCA: Sample Clustering"
- **Heatmap**: e.g. "Expression Heatmap (top N proteins)"

Use Datagrok viewer title properties or add text elements above each plot. Titles should reflect the data context (group names, number of proteins, etc.) where possible.
