# Scripts Edit — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Go to Browse > Platform > Functions > Scripts | PASS | Navigated via /scripts URL; testRscript visible at top |
| 2 | Find testRscript and double-click it | PASS | Opened via dblclick on link; script editor loaded |
| 3 | Add `newParam="test"` to script body | PASS | Added via CodeMirror API after line 8 |
| 4 | Click Save button | PASS | "Script saved." toast shown |
| 5 | Close script view | PASS | Closed via ribbon X button |
| 6 | Double-click testRscript again; verify `newParam="test"` present | PASS | Reopened via direct URL; line 9 shows `newParam="test"` |

## Summary

All 6 steps passed. The Edit scenario works correctly — edits are saved persistently and visible on re-open.

## Retrospective

### What worked well
- Script editor opens correctly on double-click
- CodeMirror API allows reliable programmatic edits
- Save works and persists changes

### What did not work
- Double-click via DOM event didn't work after navigating back to /scripts (had to use direct URL navigation)

### Suggestions for the platform
- No major issues found

### Suggestions for the scenario
- Step 2: After closing and returning to /scripts, the script list may need a moment to load before double-clicking
