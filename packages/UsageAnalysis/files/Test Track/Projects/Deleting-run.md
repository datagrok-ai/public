# Deleting — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Find projects from previous steps | PASS | PASSED | All 4 test projects found in Browse > Dashboards |
| 2 | Right-click and select Delete project | AMBIGUOUS | SKIPPED | Context menu via right-click did not appear; used grok.dapi.projects.delete() API instead |
| 3 | Click DELETE in confirmation dialog | SKIP | SKIPPED | Used API deletion, no dialog appeared |
| 4 | Check project is deleted | PASS | PASSED | All 4 test projects deleted, no longer returned by API queries |

## Summary

2 of 4 steps passed, 1 skipped, 1 ambiguous. Project deletion works via the API (`grok.dapi.projects.delete()`). The right-click context menu with "Delete project" option could not be triggered via browser automation. All test projects were successfully cleaned up.

## Retrospective

### What worked well
- API-based project deletion works reliably
- Project cleanup is idempotent (can be re-run safely)

### What did not work
- Right-click context menu on gallery items does not respond to dispatched MouseEvent('contextmenu')
- DELETE confirmation dialog was not tested (API bypass)

### Suggestions for the platform
- Ensure right-click context menus are accessible via standard DOM events
- Consider adding keyboard shortcut (e.g., Delete key) for selected project deletion

### Suggestions for the scenario
- Note says "make sure this is the last test case" — add an explicit dependency order
- Include verification step: refresh the page and confirm project is gone
