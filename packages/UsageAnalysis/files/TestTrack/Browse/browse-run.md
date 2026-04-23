# Browse — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open Browse | PASS | Browse tree opened with all sections: My stuff, Spaces, Apps, Files, Dashboards, Databases, Platform |
| 2 | Check tree items: My stuff, Spaces, Apps, Files, etc. | PASS | All sections verified. My stuff (9 items), Spaces (6 spaces), Apps (31 apps), Files (4 groups, Demo has 46 items), Databases (26 connections), Platform (5 items) |
| 3 | Run Demos | PASS | Demo app opened at /apps/Tutorials/Demo. Cheminformatics > Overview loaded: 10,000 rows x 30 cols with molecule structures. Zero console errors |
| 4 | Copy URL and open in new tab | PASS | Each section has distinct URL. Double-clicking beer.csv opened table (33 cols, 118 rows). Navigating directly to URL restored same table |
| 5 | Browse refresh remembers preview | PARTIAL | Folder-level navigation persists on refresh. But selected item preview (f=chem parameter) is dropped on reload — right panel reverts to connection details |

## Summary

4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. Item-level selection within a folder is NOT preserved across page refresh — the `?f=` query parameter is dropped on reload.

## Retrospective

### What worked well
- Browse tree shows all expected sections with correct content
- Demo app loads and runs without errors
- URL-based navigation correctly restores files and sections
- Folder-level browse state persists across refresh

### What did not work
- Item selection within a folder (the `?f=` parameter) is not preserved on page refresh — right panel reverts to parent connection details

### Suggestions for the platform
- Preserve the `?f=` query parameter across page refresh so item-level selection is maintained
- Consider adding the selected item to browser history state

### Suggestions for the scenario
- Step 1 "Open Browse" is implicit and not separately testable — consider merging with step 2
- Add explicit sub-steps for verifying expected app count and database connection count
