# Info Panels — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open linked datasets (peptides.csv) | PASS | PASSED | Loaded via `grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv')`; 647 rows, 3 cols (ID, AlignedSequence, IC50) |
| 2 | Make sure amino acids rendered with different colors | PASS | PASSED | AlignedSequence column displays colored amino acid letters (A=green, N=orange, T=teal, Y=red, K=purple, etc.) in the grid |
| 3 | Click the peptides column title | PASS | PASSED | Used `grok.shell.o = col` after enabling context panel (`grok.shell.windows.showContextPanel = true`; was hidden by default showing Help instead) |
| 4 | Check necessary panels on Context Panel (Details, Peptides) | PASS | PASSED | Context Panel shows: Details, Filter, Actions, Colors, Style, Settings, Plots, Advanced, Dev, **Peptides**, Bioinformatics |
| 5 | Expand each tab | PASS | PASSED | Expanded Details and Peptides panes via clicking pane headers |
| 6 | Make sure content for each panel displays correctly | PASS | PASSED | Details: 464 categories, weblogo chart, metadata (units, separator, aligned=SEQ.MSA, alphabet=UN, semantic type=Macromolecule); Peptides: weblogo chart positions 1-11, Activity=IC50, Scaling=none, Clusters dropdown, Generate clusters, LAUNCH SAR button |

## Summary

All 6 steps passed. The peptides.csv dataset opened correctly with amino acids color-coded in the AlignedSequence column. The Context Panel showed the expected Details and Peptides panels with correct content. Key finding: the context panel was hidden by default (showing Help instead) — enabling it via `grok.shell.windows.showContextPanel = true` was required.

## Retrospective

### What worked well
- `grok.dapi.files.readCsv()` reliably loads the demo dataset
- Amino acid color rendering is visually correct in the grid
- `grok.shell.o = col` triggers context panel content update once the panel is visible
- Both Details and Peptides panels have rich, correct content

### What did not work
- Canvas-based column header click via `MouseEvent` dispatch didn't work — the grid doesn't respond to synthetic mouse events on the canvas
- Context panel was hidden (showing Help panel instead) — not obvious from the UI

### Suggestions for the platform
- `grok.shell.windows.showContextPanel` should default to `true` or provide a more obvious way to enable it from the UI
- Add `data-column-name` attributes to column header elements for easier automation targeting

### Suggestions for the scenario
- Step 3 says "Click the peptides column title" — should clarify it's the AlignedSequence column header
- Step 4 should mention enabling the Context Panel if it's not visible (gear icon → Show Context Panel)
- Step 6 could list specific expected fields in Details (semantic type = Macromolecule) and Peptides (Activity, Scaling, Clusters, Launch SAR)
