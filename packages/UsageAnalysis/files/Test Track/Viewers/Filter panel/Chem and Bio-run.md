# Chem and Bio — Run Results

**Date**: 2026-04-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

### Section 1: Chem filter

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open spgi-100 | PASS | PASSED | 100 rows, 88 columns. Waited for Chem cell rendering (canvas poll + 3s) |
| 2 | Open the Filter Panel | PASS | PASSED | 42 filter cards shown |
| 3 | Draw CCC(N(C)C)=O in Structure filter | PASS | PASSED | Sketcher dialog, typed SMILES, clicked OK via `[name="button-OK"]`, 14 rows |
| 4 | Click gear icon, reveal search type dropdown | PASS | PASSED | `.chem-search-options-icon` click revealed `<select>` |
| 5a | Switch to Included in | PASS | PASSED | 0 rows |
| 5b | Switch to Exact | PASS | PASSED | 0 rows |
| 5c | Switch to Similar | PASS | PASSED | 0 rows |
| 5d | Switch to Not contains | PASS | PASSED | 86 rows |
| 5e | Switch to Not included in | PASS | PASSED | 100 rows |
| 5f | Switch to Contains | PASS | PASSED | 14 rows |
| 6 | Close All | PASS | PASSED | |

### Section 2: Bio

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open peptides.csv | PASS | PASSED | 647 rows, 3 columns. semType: Macromolecule |
| 2 | Wait for AlignedSequence to render | PASS | PASSED | Canvas poll + 3s wait for Bio package init |
| 3 | Open the Filter Panel | PASS | PASSED | First open: 0 filters. Close/reopen: 3 filters (AlignedSequence, ID, IC50) |
| 4 | Enter T-T-Y-K-N-Y-V in Substructure | PASS | PASSED | 28 rows filtered |
| 5 | Close All | PASS | PASSED | |

## Summary

All steps pass. All 6 chem search type modes produce correct row counts matching scenario expectations exactly. Bio substructure filter correctly filters peptides with the T-T-Y-K-N-Y-V motif (28 rows).

## Retrospective

### What worked well
- Chem: sketcher dialog, gear icon, search type switching all work reliably
- Chem: canvas poll wait ensures molecule cells render before opening filter panel
- Bio: Substructure input filters correctly via DOM events
- OK button selector `[name="button-OK"]` avoids strict mode violation (CANCEL also has `ui-btn-ok` class)

### What did not work
- Bio: first `getFiltersGroup()` call produces 0 filter cards even after canvas wait + 3s delay. Close/reopen forces fresh `look.auto()` which picks up Bio Substructure filter. Root cause: Bio package registers its filter type asynchronously after cell rendering completes, but `getFiltersGroup()` calls `look.auto()` only on initial open

### Suggestions for the platform
- `getFiltersGroup()` should re-run `look.auto()` if no filter cards were generated, or listen for late filter type registrations from packages
- Consider emitting an event when all packages finish registering their filter types

### Suggestions for the scenario
- Bio section step 2 is critical — filter panel must be opened only after Bio package fully initializes
- The scenario doesn't mention the close/reopen workaround needed for Bio datasets
