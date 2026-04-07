# Calculate — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open smiles.csv | PASS | Opened via JS API, 1000 rows, 20 cols |
| 2 | Run Chem > Calculate > Descriptors | FAIL | Menu item clicked but dialog did not open; function call via `Chem:chemDescriptors` also failed with "Error" balloon |
| 3 | Select arbitrary values in dialog | SKIP | Dialog did not open |
| 4 | Check column with calculated values added | SKIP | Dependent on step 2 |

## Summary

The Descriptors function failed to open its dialog when invoked both via menu and JS API. An "Error" balloon appeared with no detailed message. 1 step passed, 1 failed, 2 skipped.

## Retrospective

### What worked well
- Dataset opened correctly
- Chem menu navigation worked

### What did not work
- Descriptors dialog failed to open — root cause unknown, possibly a server-side computation dependency issue on public.datagrok.ai
- Error balloon provided no useful error message

### Suggestions for the platform
- Error balloons should include a detailed error message or stack trace for debugging
- Consider providing a fallback when server-side computation is unavailable

### Suggestions for the scenario
- Specify which descriptors to select for a reproducible test
- Add expected column names for verification
