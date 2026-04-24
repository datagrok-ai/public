# Calculate — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open smiles.csv dataset | 9s | PASS | PASSED | 1000 rows, canonical_smiles Molecule col |
| 2 | Open Chem → Calculate → Descriptors... dialog | 15s | FAIL | FAILED | Top menu click emits no `.d4-menu-popup`; `grok.shell.topMenu.find('Chem').find('Calculate').find('Descriptors...').click()` no-ops too |
| 3 | Accept defaults → new descriptor columns appended | n/a | SKIP | FAILED | Could not reach OK — dialog never opened |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 6m |
| grok-browser execution (scenario steps) | 2m |
| Execute via grok-browser (total) | 8m |
| Spec file generation | 1m |
| Spec script execution | 38s (FAILED) |
| **Total scenario run (with model)** | ~10m |

## Summary

Calculate Descriptors cannot be exercised on `dev` right now. The Chem top menu fails to open its popup — both through DOM `dispatchEvent(new MouseEvent('click'))` and via `grok.shell.topMenu.find(...).click()`. Console logs show repeated `Connection to db.datagrok.ai:54325 refused` errors, which correlate with DB-backed initialization failing. R-Groups Analysis (scenario 2) succeeded earlier in the same session — suggesting the regression is specific to the Chem Calculate/Descriptors function registration or its menu publication, not the whole Chem package.

## Retrospective

### What worked well
- Dataset opening and the browse tree work fine
- Other Chem sub-commands (R-Groups Analysis) continue to function in the same environment

### What did not work
- Top menu popup never renders when `[name="div-Chem"]` is clicked via synthesized events
- `grok.shell.topMenu.find('Chem').find('Calculate').find('Descriptors...').click()` returns silently
- Dev DB cannot be reached from the server (`db.datagrok.ai:54325` refused); several package bootstraps log `Promise rejected` in console

### Suggestions for the platform
- Investigate Chem Descriptors menu registration — does the registered function silently fail its precondition (e.g., "chemDescriptorsTree" server call) when DB is unreachable?
- Menu items should still open their dialogs even when an async prerequisite (tree fetch, DB handshake) fails; degrade to showing an error balloon rather than a no-op click
- The `Connection to db.datagrok.ai:54325 refused` on dev needs an infra check — several warnings about "Failed to load DB schema, Object handlers not registered" show the same root cause

### Suggestions for the scenario
- Specify the expected dialog contents: tree of descriptor categories, a column selector for the molecule column, and an OK button
- For the "Repeat for each section under Calculate" bullet, enumerate the sections explicitly (Chemical Properties, Toxicity Risks, MPO Score, Biochemical Properties, Descriptors, Map Identifiers, …) so QA can tick them off
