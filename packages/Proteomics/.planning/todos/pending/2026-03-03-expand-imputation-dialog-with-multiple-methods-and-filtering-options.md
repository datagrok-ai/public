---
created: 2026-03-03T19:06:00.000Z
title: Expand imputation dialog with multiple methods and filtering options
area: analysis
files:
  - packages/Proteomics/src/analysis/imputation.ts
  - packages/Proteomics/src/package.ts
---

## Problem

The current imputation applies MinProb with hardcoded parameters and no user control. Users need to choose imputation methods, tune parameters, and filter proteins by minimum valid values before imputing.

## Solution

Replace the current auto-impute with a rich dialog:

```
┌─────────────────────────────────────────┐
│  Impute Missing Values                  │
├─────────────────────────────────────────┤
│  Table:        [current table    ▼]     │
│  Intensity columns: [select...   ▼]     │
│  Method:       [MinProb          ▼]     │
│                                         │
│  MinProb settings:                      │
│    Quantile:     [0.01          ]       │
│    Downshift:    [1.8  σ        ]       │
│                                         │
│  Filter proteins:                       │
│    Min valid values: [2         ]       │
│    In at least one group: [☑]           │
│    Group column: [Condition      ▼]     │
│                                         │
│            [Cancel]  [OK]               │
└─────────────────────────────────────────┘
```

Key features:
- **Method dropdown**: MinProb, KNN, zero, mean/median, MNAR-aware methods
- **Method-specific settings**: Show/hide parameter controls based on selected method (e.g. quantile + downshift for MinProb, k for KNN)
- **Protein filtering**: Min valid values threshold with option to require them within at least one experimental group
- **Group-aware filtering**: Use the annotated experiment groups to determine per-group valid counts
- **Column selection**: Let user pick which intensity columns to impute
