---
created: 2026-03-03T19:04:00.000Z
title: Expand normalize dialog with method selection and options
area: analysis
files:
  - packages/Proteomics/src/analysis/normalization.ts
  - packages/Proteomics/src/package.ts
---

## Problem

The current normalization function applies median centering automatically with no user control. Users need to choose normalization methods and configure options like log2 transform, before/after visualization, and whether to replace original columns.

## Solution

Replace the current auto-normalize with a rich dialog:

```
┌─────────────────────────────────────────┐
│  Normalize Proteomics Data              │
├─────────────────────────────────────────┤
│  Table:        [current table    ▼]     │
│  Intensity columns: [select...   ▼]     │
│  Method:       [Median centering ▼]     │
│                                         │
│  ☑ Log2 transform first                 │
│  ☑ Show before/after plots              │
│  ☐ Replace original columns             │
│                                         │
│            [Cancel]  [OK]               │
└─────────────────────────────────────────┘
```

Key features:
- **Method dropdown**: Median centering, quantile normalization, VSN, etc.
- **Log2 transform checkbox**: Option to log2 transform before normalizing
- **Before/after plots**: Show box plots or density plots comparing distributions pre/post normalization
- **Replace original**: Toggle between creating new columns vs overwriting
- **Column selection**: Let user pick which intensity columns to normalize
