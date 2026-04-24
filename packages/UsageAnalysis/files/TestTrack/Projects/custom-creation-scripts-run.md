# Custom creation scripts — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: SKIP

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Run the getLastCreatedFile script | SKIP | SKIPPED | Script execution and Data Sync not tested |
| 2 | Add viewers and save with Data Sync | SKIP | SKIPPED | Dependent on step 1 |
| 3 | Close all | SKIP | SKIPPED | Dependent on step 1 |
| 4 | Update CSV file in DemoFiles/chem | SKIP | SKIPPED | Dependent on step 1 |
| 5 | Open saved project, verify data refreshed | SKIP | SKIPPED | Dependent on step 1 |

## Summary

All 5 steps skipped. This scenario requires running a custom JavaScript script with Data Sync enabled, then modifying files on the server to verify that the project reloads fresh data. Not feasible without dedicated server-side file manipulation capabilities.

## Retrospective

### What worked well
- N/A (scenario not executed)

### What did not work
- Script execution with Data Sync requires server-side setup

### Suggestions for the platform
- Provide test utilities for script-based Data Sync verification

### Suggestions for the scenario
- Add pre-condition: which file to create/modify in chem folder
- Clarify expected behavior when no numeric suffix files exist
