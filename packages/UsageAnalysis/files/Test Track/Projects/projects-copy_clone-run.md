# Projects copy_clone — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: SKIP

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Browse and check preview of all created projects | PASS | PASSED | All test projects visible with thumbnails in gallery view |
| 2 | Share created projects | SKIP | SKIPPED | Context menu not accessible via automation |
| 3 | Open created projects | PASS | PASSED | Projects open correctly via API |
| 4 | Edit, save copy with link, save copy with clone | SKIP | SKIPPED | Save dialog copy/link/clone options require complex UI interaction |
| 5 | Share new projects, check ability to open shared | SKIP | SKIPPED | Dependent on step 4 |

## Summary

2 of 5 steps passed, 3 skipped. Project preview and opening work. Copy/clone/link operations were not tested because the Save dialog's advanced options (Save a copy, Clone, Link) require multi-step UI interactions that were not feasible via current automation.

## Retrospective

### What worked well
- Project gallery previews render correctly
- Projects open and display data correctly

### What did not work
- Copy/Clone/Link dialog options were not accessible for automation
- Sharing via context menu could not be triggered

### Suggestions for the platform
- Expose project copy/clone operations via JS API for automation
- Add `data-testid` to Save dialog options (Clone/Link/Move radio buttons)

### Suggestions for the scenario
- Break down step 4 into sub-steps with expected UI state after each action
- Specify exact viewer to add (e.g., "add a Scatter plot")
