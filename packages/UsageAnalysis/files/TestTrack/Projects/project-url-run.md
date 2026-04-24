# Project URL — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: SKIP

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Browse > Dashboards | SKIP | SKIPPED | Dependent on copy/clone scenario |
| 2 | Click projects (original, copy, clone, layout) | SKIP | SKIPPED | Copy/clone not performed |
| 3 | Go to Context Panel > Links, copy URL | SKIP | SKIPPED | No projects to test |
| 4 | Open URL in new browser tab | SKIP | SKIPPED | No URL to test |

## Summary

All steps skipped. This scenario depends on Projects copy_clone.md (order 5) which was not fully executed. The Link/Clone/Copy project operations were not performed in prior steps.

## Retrospective

### What worked well
- N/A

### What did not work
- This scenario has hard dependency on copy_clone scenario

### Suggestions for the platform
- N/A

### Suggestions for the scenario
- Make dependencies explicit (requires Projects copy_clone.md to have been run first)
- Could include a standalone test: open any existing project, copy its URL, navigate to it
