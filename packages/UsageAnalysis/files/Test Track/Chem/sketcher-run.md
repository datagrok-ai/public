# Sketcher — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open smiles.csv and show sketcher dialog with first molecule | 6s | PASS | PASSED | Sketcher loads via `grok.chem.sketcher(null, smiles)` + `ui.dialog(...)` |
| 2 | Enter `C1CCCCC1` into molecular SMILES input + Enter | 3s | PASS | PASSED | Input accepts value and fires onInput |
| 3 | Hamburger menu exposes Copy as SMILES / MOLBLOCK / Recent / Favorites | 2s | PASS | PASSED | Expected menu items present in `.d4-menu-item-label` list |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 30s |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | 30s |
| Spec file generation | 30s |
| Spec script execution | 16.5s |
| **Total scenario run (with model)** | ~1m 30s |

## Summary

Sketcher opens via `grok.chem.sketcher(molCol, initialSmiles)` wrapped in `ui.dialog(...).show()`, accepts a typed SMILES (`C1CCCCC1`), and the hamburger menu exposes the expected Copy as SMILES / Copy as MOLBLOCK / Recent / Favorites options. Deep interactions (double-click on a grid cell, Favorites add/remove, per-sketcher-type iteration) remain manual because sketcher canvas is not DOM-inspectable.

## Retrospective

### What worked well
- `grok.chem.sketcher(null, smiles)` returns a widget that plugs into `ui.dialog()` cleanly
- SMILES text input works via native value setter + `input`/`keydown` dispatch

### What did not work
- Scenario's "double-click a molecule in the grid" can't be reproduced in automation — the grid cell is on canvas; invoking the sketcher directly via API is the only viable replay
- Iterating through every sketcher type (Marvin, Chem Draw, OpenChemLib, Ketcher) requires the hamburger → Sketcher menu, which is user-specific and not scriptable without deep setup

### Suggestions for the platform
- Expose a sketcher-switcher API like `grok.chem.setDefaultSketcher('ketcher')` so automation can iterate sketcher types
- Persist the last-used sketcher type per user setting (not per session) — automation currently has no stable handle on which sketcher is active

### Suggestions for the scenario
- Move the "check all sketcher types" bullet to its own sub-scenario; it's a matrix test
- `CC[C@H]...` giant peptide bug reference (#1608) deserves its own scenario so regression tests have a single clear assertion
