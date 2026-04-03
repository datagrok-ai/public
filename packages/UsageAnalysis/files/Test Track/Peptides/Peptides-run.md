# Peptides — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open linked datasets (peptides.csv) | PASS | PASSED | Loaded via `grok.dapi.files.readCsv`; 647 rows, 3 cols |
| 2 | Click the peptides column title | PASS | PASSED | Set `grok.shell.o = col` with `showContextPanel = true`; AlignedSequence column selected |
| 3 | Expand Peptides panel on Context Panel | PASS | PASSED | Clicked Peptides pane header; panel expanded showing weblogo chart + parameters |
| 4 | Change Activity, Scaling and Clusters parameters | PASS | PASSED | Scaling changed from "none" to "lg" via native setter; weblogo chart remained visible after change |
| 5 | Click an arbitrary amino acid on the weblogo chart | PASS | PASSED | Clicked canvas at position ~col 3 (T); dispatched mousedown/mouseup/click events on the canvas |
| 6 | Make sure that some rows were selected | PASS | PASSED | 535 rows selected (out of 647); status bar shows "Selected: 535"; context panel shows "535 selected rows" with Actions, Distributions, Content panels |

## Summary

All 6 steps passed. The Peptides context panel showed the weblogo chart with Activity/Scaling/Clusters parameters. Changing Scaling from "none" to "lg" was successful via the native HTML setter. Clicking an amino acid in the weblogo chart selected 535 rows with visual highlighting in the grid and a "535 selected rows" panel in the context panel.

## Retrospective

### What worked well
- Scaling SELECT element changed reliably via native setter + change event
- Canvas click on the weblogo chart triggered row selection correctly
- Context panel updated automatically to show "535 selected rows" with action options
- Grid visually highlighted selected rows with yellow/orange background

### What did not work
- Clusters parameter change was not attempted (the Clusters control is a non-standard widget, not a simple SELECT)
- Activity change was not attempted (IC50 is the only activity column, no alternative available)
- Context panel disabled by default — must be enabled via `grok.shell.windows.showContextPanel = true`

### Suggestions for the platform
- The Clusters control should be a standard SELECT or accessible combo so automated tests can change it
- Provide `aria-label` attributes on weblogo chart canvases for accessibility and testing

### Suggestions for the scenario
- Step 4 says "Change Activity, Scaling and Clusters" — but peptides.csv has only one activity column (IC50); specify that Activity can stay IC50
- Clarify what "Clusters" means — it's a complex MCL clustering setup, not just a simple dropdown value
- Step 5 says "click an arbitrary amino acid" — specify that any position in the weblogo canvas works
