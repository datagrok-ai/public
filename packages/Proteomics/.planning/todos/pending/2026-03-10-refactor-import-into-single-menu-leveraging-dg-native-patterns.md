---
created: 2026-03-10T14:16:05.170Z
title: Refactor import into single menu leveraging DG native patterns
area: import
files:
  - packages/Proteomics/src/package.ts
  - packages/Proteomics/src/package.g.ts
---

## Problem

The Proteomics package currently registers separate menu items for each import format (MaxQuant, Spectronaut, Generic Matrix). This creates a fragmented import experience. Datagrok has native patterns for unified file handling -- file handlers, auto-detection, and file info panels -- that would provide a more natural UX where the platform routes files to the correct parser automatically based on content sniffing or extension.

## Solution

Refactor to a single "Import Proteomics Data" entry point (or leverage Datagrok's `//meta.ext` and `//input: file` patterns) that:
1. Auto-detects format from file content (MaxQuant proteinGroups.txt headers, Spectronaut required columns, generic CSV/TSV fallback)
2. Uses Datagrok's native file handler registration so files are recognized on drag-and-drop or open
3. Falls back to a format selector dialog if auto-detection is ambiguous
4. Consolidates the three separate import handlers into a single dispatch function

See also: existing todo for expanding import sources (local drive, DG folders, existing DataFrames) -- these should be combined in the refactor.
