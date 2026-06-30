---
created: 2026-03-03T19:08:00.000Z
title: Expand DE dialog with method selection, comparison picker, and options
area: analysis
files:
  - packages/Proteomics/src/analysis/differential-expression.ts
  - packages/Proteomics/src/package.ts
---

## Problem

The current DE dialog has basic group info and threshold inputs. Users need richer controls: method selection, dynamic comparison from group annotations, optional volcano plot auto-open, and peptide count column picker for DEqMS.

## Solution

Replace the current DE dialog with a full-featured version:

```
┌─────────────────────────────────────────────┐
│  Differential Expression Analysis           │
├─────────────────────────────────────────────┤
│  Table:           [current table      ▼]    │
│  Intensity cols:  [select...          ▼]    │
│  Group column:    [Condition           ▼]   │
│                                             │
│  Comparison:      [Treatment vs Control ▼]  │
│    (auto-populated from group column)       │
│                                             │
│  Method:          [DEqMS               ▼]   │
│  Peptide count col: [Peptide counts    ▼]   │
│    (only shown when DEqMS selected)         │
│                                             │
│  FDR threshold:   [0.05               ]     │
│  Log2FC threshold: [1.0               ]     │
│                                             │
│  ☑ Add volcano plot                         │
│  ☑ Add results to table                     │
│                                             │
│              [Cancel]  [OK]                 │
└─────────────────────────────────────────────┘
```

Key features:
- **Group column picker**: Auto-detect from experiment annotation, populate comparison dropdown from unique values
- **Comparison dropdown**: Auto-generate pairwise comparisons from group levels (Treatment vs Control, etc.)
- **Method dropdown**: limma, DEqMS, t-test (client) — with conditional fields per method
- **Peptide count col**: Only visible when DEqMS selected (already partially implemented)
- **Volcano plot checkbox**: Auto-open volcano viewer after DE completes
- **Add results checkbox**: Toggle whether log2FC/p-value columns are added to the table
- Subsumes the existing "add local t-test" todo — t-test becomes one of the method options

## Notes
- This supersedes the narrower "add local t-test back" todo once implemented
