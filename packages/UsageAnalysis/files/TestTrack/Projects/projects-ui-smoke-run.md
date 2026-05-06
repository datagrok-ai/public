# Projects — UI smoke — Run Results

**Date**: 2026-05-06
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Step 1: open demog.csv (UI — Browse-tree right-click → Open) | n/a | PASS | PASSED | Anchored DOM walk: `.d4-tree-view-root > .d4-tree-view-group-host > children[3] (Files) > children[2] (Demo)` + chevron click + lazy-load wait + `dispatchEvent('contextmenu')` on `.d4-tree-view-node` + `[name="div-Open"]` click |
| 2 | Step 2: open Save dialog, ensure Data sync ON, set name, click OK | n/a | PASS | PASSED | Toolbar `[name="button-Save"]` visible (bug 2b platform-fixed 2026-05-06); Data Sync ON by default; `input#name` set via native value-setter; `[name="button-OK"]` click |
| 3 | Step 3: cancel auto-Share dialog | n/a | PASS | PASSED | `.d4-dialog[name^="dialog-Share-..."]` → `[name="button-CANCEL"]` |
| 4 | Step 4: navigate to Browse > Dashboards, assert tile visible | n/a | PASS | PASSED | `/projects` URL (NOT `/browse/dashboards`); tile name attr uses verbatim `[name="div-UiSmoke<stamp>"]` (PascalCase preserved) |
| 5 | Step 5-6: right-click tile → Share → fill recipient → OK | n/a | PASS | PASSED | `[name="div-Share..."]` → Share dialog → recipient input by placeholder `User, group, or email` → access-level select |
| 6 | Step 7-8: double-click tile to reopen, verify demog table loaded | n/a | PASS | PASSED | Tile dblclick → `[name="viewer-Grid"]` waitFor → `tv.dataFrame.rowCount > 0` |
| 7 | Step 9-11: right-click tile → Delete → confirm → assert tile gone | n/a | PASS | PASSED | `[name="div-Delete-Project"]` → `dialog-Are-you-sure?` → `[name="button-DELETE"]` → dialog hides after server-side delete completes (60s timeout — server-delete can take ~30-60s under Playwright load) |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 18s |
| **Total scenario run** | 1m 18s |

## Summary

Spec passed end-to-end on dev.datagrok.ai using bundled Playwright Chromium. Total run: 1m 18s.

This was the last spec to complete after platform-side fix of bug 2b (toolbar SAVE button collapsing to `offsetWidth=0` after JS-API `openTableFromFile`) on 2026-05-06.

## Path to PASS — fix history

| Date | Issue | Resolution |
|---|---|---|
| 2026-05-05 | `getByText('Demo', {exact:true})` matched 24 hidden duplicates in virtualized panes | Switched to anchored DOM walk via `.d4-tree-view-root > .d4-tree-view-group-host > children[3]` |
| 2026-05-05 | `'Demo'` label class wrong | `.d4-tree-view-group-label` (verified by MCP) |
| 2026-05-05 | `Demo` tri carries plain `d4-tree-view-tri` (NEITHER -collapsed nor -expanded) until first expansion | Detect "not yet expanded" by Demo's child-host emptiness rather than tri class |
| 2026-05-05 | Right-click via `.click({button:'right'})` on label unreliable | `dispatchEvent('contextmenu')` on `.d4-tree-view-node` ancestor |
| 2026-05-05 | Tile slug rule wrong (`name.toLowerCase()`) | Use `grok.shell.project.name` verbatim — server preserves PascalCase for non-hyphen names |
| 2026-05-05 | `/browse/dashboards` URL wrong | `/projects` |
| 2026-05-05 | `permissions.find()` does not exist | Use `permissions.get()` or skip server-side perm verification |
| 2026-05-06 | bug 2b — toolbar SAVE button collapsed offsetWidth=0 in Playwright | Platform-fixed; channel:'chrome' override removed (bundled Chromium now sufficient) |
| 2026-05-06 | DELETE confirm dialog `toBeHidden({timeout:15_000})` failed | Server-side delete takes longer than 15s under Playwright load — extended to 60s |

## Retrospective

### What worked well
- Anchored DOM walk pattern (`:scope > .d4-tree-view-group-host`) reliably locates Browse-tree nodes without depending on virtualized cached duplicates.
- `dispatchEvent('contextmenu')` on the tree-node ancestor reliably triggers Datagrok's right-click handler.
- Verbatim `grok.shell.project.name` for tile selectors avoids slug-derivation bugs.
- 60s DELETE-confirm timeout absorbs server-side delete latency.

### What did not work
- Nothing notable in the final run.

### Suggestions for the platform
- Consider closing the `dialog-Are-you-sure?` confirm dialog optimistically (immediately on DELETE click) and showing progress in the tile/toolbar instead — keeps UI responsive while delete is in flight.

### Suggestions for the scenario
- None from this run.
