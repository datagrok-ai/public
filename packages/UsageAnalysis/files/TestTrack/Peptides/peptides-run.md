# Peptides — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open peptides.csv | PASS | 8s | PASSED | 647 rows loaded, AlignedSequence detected as Macromolecule |
| 2 | Click peptides column title | PASS | 2s | PASSED | grok.shell.o = column used to trigger Context Panel update |
| 3 | Expand Peptides panel on Context Panel | PASS | 2s | PASSED | Peptides pane shows Activity (IC50), Scaling, Clusters, Launch SAR button, WebLogo and Histogram |
| 4 | Change Activity, Scaling and Clusters parameters | PASS | 3s | PASSED | Scaling changed to 'lg', Generate clusters checkbox toggled |
| 5 | Click amino acid on weblogo chart | AMBIGUOUS | - | SKIP | WebLogo canvas exists (2 canvases in Peptides pane) but dispatched click events don't trigger amino acid selection |
| 6 | Verify some rows were selected | AMBIGUOUS | - | SKIP | Depends on step 5; selection count remained 0 after dispatched clicks |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 18s |
| Spec file generation | 3s |
| Spec script execution | 11.6s |

## Summary

Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clusters parameters and a WebLogo chart. Parameter changes work via DOM events. Steps 5-6 are ambiguous because the WebLogo canvas does not respond to dispatched mouse events for amino acid selection.

## Retrospective

### What worked well
- Peptides.csv loads with Macromolecule semType detection
- Peptides pane in Context Panel shows all expected controls (Activity, Scaling, Clusters, Launch SAR)
- WebLogo canvas and Histogram canvas are present in the Peptides pane
- Parameter changes (Scaling dropdown, Generate clusters checkbox) work via DOM events
- `grok.shell.o = column` reliably triggers Context Panel update for column-level panels

### What did not work
- WebLogo canvas click does not trigger amino acid selection via dispatched DOM events
- Context Panel loses column focus when canvas is scrolled into view — reverts to table-level panes
- `df.currentCol = col` does not reliably trigger column-level Context Panel display

### Suggestions for the platform
- Expose WebLogo amino acid selection API for programmatic use
- Prevent Context Panel from losing focus when scrolling within a pane

### Suggestions for the scenario
- Note that WebLogo interaction requires real browser clicks (cannot be automated via JS events)
- Specify which amino acid position to click for reproducibility
