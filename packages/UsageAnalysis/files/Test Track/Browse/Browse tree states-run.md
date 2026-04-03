# Browse Tree States — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Browse tree should remember the last expand/collapse state | PARTIAL | **Within session**: PASS — expanded Spaces, opened table view via grok.data.demo.demog(), closed all, returned to Browse: Spaces still expanded with 7 children visible. **Across page refresh**: FAIL — after F5 reload, Spaces collapsed back to default state. Tree remembers state within a session but not across page refreshes. |

## Summary

1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a single session (when switching between table views and Browse). However, the tree state is NOT persisted across page refreshes — nodes revert to a default expanded state after reload.

## Retrospective

### What worked well
- Tree state is preserved within a session when switching between views
- Default expand state on load is reasonable (My stuff, Files, Databases, Platform expanded)

### What did not work
- Expand/collapse state is not persisted across page refreshes — expanded "Spaces" node collapsed after reload

### Suggestions for the platform
- Persist tree expand/collapse state in localStorage or server-side user preferences so it survives page refresh
- Consider saving scroll position in the tree as well

### Suggestions for the scenario
- Clarify what "last expand/collapse state" means — within a session? Across page refresh? Across browser restart?
- Add explicit sub-steps: (a) expand a collapsed node, (b) switch views, (c) return to Browse and verify, (d) refresh page and verify
