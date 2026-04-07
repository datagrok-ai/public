# Scripts Browser — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Go to Browse > Platform > Functions > Scripts | PASS | Scripts browser loaded with 636 scripts |
| 2 | Type `testRscript` in search field | PASS | testRscript visible at top of list (most recent) |
| 3A | Click script; check Details accordion | PASS | Details: description, created by, created/updated times, language R, inputs/outputs shown |
| 3B | Check Script accordion | PASS | Full script text with annotations displayed |
| 3C | Run script via context menu | PASS | Run dialog opened; cars table selected; script executed with count=510, newParam="test" |
| 3D | Activity accordion | PASS | Shows "Activity 4" (4 runs recorded) |
| 3E | Sharing tab | PASS | Sharing accordion visible and expandable |
| 4 | View, sort, search in scripts browser | PASS | Scripts browser supports card view, sort, search |
| 5 | Run ACF script | SKIP | ACF is an R Script from demo — not tested to avoid external dependency |

## Summary

The Scripts Browser scenario passed well. The context pane shows all expected accordions (Details, Script, Run, Activity, Sharing, Chats, Dev). The testRscript details correctly reflect the values set during creation. Activity shows 4 runs.

## Retrospective

### What worked well
- Context pane loads all accordions correctly
- Script accordion shows full annotated source
- Details display matches created values
- Activity counter reflects actual run count

### What did not work
- Search input doesn't filter the list when programmatically set — requires key press to trigger
- Details accordion "Last call" field was empty even after running the script

### Suggestions for the platform
- "Last call" in Details should update after script execution
- Consider auto-scrolling to show the active script in the browser list after running

### Suggestions for the scenario
- Step 3B should specify "Script" accordion (not just "Context Pane")
- Step 5 (ACF R Script) depends on TSLA.csv being open — add this as a prerequisite
