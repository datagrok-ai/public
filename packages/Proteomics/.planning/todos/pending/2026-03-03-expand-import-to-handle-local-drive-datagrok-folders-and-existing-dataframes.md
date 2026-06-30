---
created: 2026-03-03T19:02:01.219Z
title: Expand import to handle local drive, Datagrok folders, and existing dataframes
area: import
files:
  - packages/Proteomics/src/parsers/maxquant-parser.ts
  - packages/Proteomics/src/package.ts
---

## Problem

The current import function only handles file upload (drag-and-drop / file open). Users should also be able to load proteomics datasets from:
1. **Local drive** — browse and select a proteinGroups.txt from disk
2. **Datagrok folders** — pick from files already uploaded to the platform (e.g. `System:AppData/Proteomics/` or user file shares)
3. **Existing dataframes** — apply the proteomics pipeline to a DataFrame already open in the platform

The screenshot shows a UI pattern with icons for these sources (dropdown, folder browse, Datagrok files, existing tables).

## Solution

Expand the dataset input to support multiple sources. Use Datagrok's `ui.input.file()` or file browser APIs for local/server files. For existing dataframes, use `ui.input.table()` or `grok.shell.tables`. Wire all sources through the same MaxQuant parser / column detection pipeline.
