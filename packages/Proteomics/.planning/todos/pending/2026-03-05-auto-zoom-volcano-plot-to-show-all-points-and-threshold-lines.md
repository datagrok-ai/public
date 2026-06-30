---
created: 2026-03-05T00:00:00.000Z
title: Auto-zoom volcano plot to show all points and threshold lines
area: ui
files:
  - packages/Proteomics/src/viewers/volcano.ts
---

## Problem

The volcano plot's default axis range may not include all data points or the significance/fold-change threshold lines. Users have to manually zoom or pan to see the full picture, which is a poor default experience.

## Solution

When creating the volcano plot, compute axis ranges that encompass:
1. All data points (min/max of log2FC and -log10 p-value)
2. The FC threshold lines (±fcThreshold on X axis)
3. The p-value threshold line (-log10(pThreshold) on Y axis)
4. Add a small padding margin (e.g. 5-10%) beyond the extremes

Set the scatter plot's X/Y min/max properties accordingly so the default view shows everything without manual adjustment.
