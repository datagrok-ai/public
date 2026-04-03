# Input Functions — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: SKIP

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open spgi.csv dataset | SKIP | N/A | spgi.csv not available as demo table |
| 2 | Open Add New Column dialog | SKIP | N/A | Depends on step 1 |
| 3 | Add function via plus icon | SKIP | N/A | Depends on step 2 |
| 4 | Add function via drag-and-drop | SKIP | N/A | Depends on step 2 |
| 5 | Add column via drag-and-drop | SKIP | N/A | Depends on step 2 |
| 6 | Verify formula text field content | SKIP | N/A | Depends on previous steps |

## Summary

All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on release-ec2.

## Retrospective

### What worked well
- N/A (scenario not executed)

### What did not work
- spgi.csv dataset not available on the test server

### Suggestions for the platform
- Include spgi.csv in the standard demo datasets

### Suggestions for the scenario
- Specify how to obtain/load spgi.csv, or use demog.csv which is universally available
- Drag-and-drop interactions are difficult to automate — consider providing keyboard/menu alternatives for testing
