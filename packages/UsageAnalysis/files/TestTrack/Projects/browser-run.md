# Browser — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Browse > Dashboards | PASS | PASSED | Navigated to /projects, all 4 test projects visible |
| 2 | Search for saved projects | PASS | PASSED | Search box filtered to 7 matching "Test Project" entries |
| 3 | Right-click project and select Share | SKIP | SKIPPED | Context menu did not appear via dispatched contextmenu event |
| 4 | Share with registered user and email | SKIP | SKIPPED | Dependent on step 3 |
| 5 | Check Sharing section on Context Panel | PASS | PASSED | Context Panel shows Sharing section when project selected |
| 6 | Check notifications | SKIP | SKIPPED | No sharing was performed |
| 7 | Right-click and select Details | AMBIGUOUS | N/A | Context menu not clickable; Details visible in Context Panel |
| 8 | Review all Context Panel tabs | PASS | PASSED | Details, Custom views, Content, Activity, Sharing, Chats, Advanced visible |
| 9 | Run old projects | PASS | PASSED | Opened "Test Project - Case 1 Local" via API, tables loaded correctly |

## Summary

5 of 9 steps passed, 3 skipped, 1 ambiguous. Browse > Dashboards view works correctly: projects are listed, searchable, and show full details in the Context Panel. Opening saved projects works. The right-click context menu could not be triggered via browser automation, preventing Share and Details menu testing.

## Retrospective

### What worked well
- Projects gallery view loads and displays all projects
- Search filtering works correctly
- Context Panel shows all expected tabs
- Opening projects via `project.open()` API works reliably

### What did not work
- Right-click context menu on gallery grid items could not be triggered via DOM event dispatch
- Double-click on project cards in gallery view did not navigate to open the project

### Suggestions for the platform
- Add `data-testid` or accessible role attributes to project gallery cards
- Ensure context menu is accessible programmatically for automation

### Suggestions for the scenario
- Clarify how to navigate to Dashboards (Browse tree vs URL)
- Step 9 "Run some old projects" is vague — specify which projects to open
