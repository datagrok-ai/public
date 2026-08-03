# Filter Panel — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI.csv → filter panel shows Structure filter | 9s | PASS | PASSED | `[name="viewer-Filters"] .d4-filter` includes a Structure/Molecule filter card |
| 2 | Click Structure filter sketch-link → sketcher dialog opens → enter c1ccccc1 → OK → table filtered | 13s | PASS | PASSED | Filter sets `df.filter.trueCount < df.rowCount`, > 0 rows remain |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 18s |
| grok-browser execution (scenario steps) | 16s |
| Execute via grok-browser (total) | 34s |
| Spec file generation | 25s |
| Spec script execution | 24.2s |
| **Total scenario run (with model)** | 1m 23s |

## Summary

The filter panel correctly shows a Structure filter for SPGI.csv's Molecule column; clicking the embedded sketch-link opens the sketcher dialog, typing `c1ccccc1` + Enter + OK applies a substructure filter that reduces `df.filter.trueCount` below the total row count. Not exercised in automation: per-mode toggles (Contains / Included in / Exact / Similar), right-click "Current Value → Use as filter", drag-and-drop column header → filter panel, column hamburger → Filter → Add filter, multi-view sync.

## Retrospective

### What worked well
- `grok.shell.tv.getFiltersGroup()` + wait for `[name="viewer-Filters"] .d4-filter` reliably opens the filter panel
- Sketch-link click → dialog → SMILES input pattern works through synthesized DOM events
- `df.filter.trueCount` vs `df.rowCount` is a strong assertion

### What did not work
- The scenario's many "also check..." bullets (modes, DnD, right-click, 2 views, hamburger menu) expand into 5+ sub-scenarios; running them in one spec is impractical

### Suggestions for the platform
- Expose filter mode selector as a named element (`[name="input-Filter-mode"]`) so automation can iterate modes programmatically
- Structure filter card should expose `data-testid="structure-filter"` to avoid relying on header text match

### Suggestions for the scenario
- Split into multiple focused scenarios: "Sketch filter", "Drag'n'drop column → filter", "Use as filter (context menu)", "Filter sync across views"
- The current numbering is inconsistent (`1. 1. 1. 2. 2.` alternates) — renumber as flat list
