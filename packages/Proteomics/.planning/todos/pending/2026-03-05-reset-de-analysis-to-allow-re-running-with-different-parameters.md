---
created: 2026-03-05T00:00:00.000Z
title: Reset DE analysis to allow re-running with different parameters
area: analysis
files:
  - packages/Proteomics/src/analysis/differential-expression.ts
  - packages/Proteomics/src/package.ts
---

## Problem

Once DE is run, the `proteomics.de_complete` tag is set and `showDEDialog()` blocks re-running with "Differential expression already performed". Users need to iterate with different thresholds, methods, or group assignments without reimporting the data.

## Solution

Add a "Reset DE" action (menu item or button in the DE dialog) that:
1. Removes the `log2FC`, `p-value`, `adj.p-value`, and `significant` columns
2. Clears the `proteomics.de_complete` and `proteomics.de_method` tags
3. Fires `valuesChanged` so viewers update
4. Allows the user to re-open the DE dialog with new parameters

Could also auto-detect existing DE columns in `showDEDialog()` and offer to replace them instead of blocking.
