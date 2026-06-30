---
created: 2026-03-03T19:38:00.000Z
title: Table should retain the name of the imported file
area: import
files:
  - packages/Proteomics/src/parsers/maxquant-parser.ts
---

## Problem

When importing a proteinGroups.txt file, the resulting DataFrame doesn't preserve the original filename. Users lose track of which file produced which table, especially when working with multiple datasets.

## Solution

Set the DataFrame name to the imported filename (e.g. "proteinGroups" or "cptac-spike-in") during parsing. Use `df.name = filename` (sans extension) in the MaxQuant parser after creating the DataFrame.
