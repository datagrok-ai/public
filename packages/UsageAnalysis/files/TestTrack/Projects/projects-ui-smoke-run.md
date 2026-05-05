# Projects — UI smoke — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: FAIL (in-progress; pending MCP-driven Step 1/2 selector verification in new session)

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Step 1: open demog.csv (UI — Browse-tree right-click → Open) | n/a | FAIL | FAILED | Anchored DOM walk + chevron click + right-click → `[name="div-Open"]` runs without throwing, but `grok.shell.tv?.dataFrame` is `undefined` afterwards — `expect(info.rows).toBeGreaterThan(0)` receives `undefined` |
| 2 | Step 2: open Save dialog, ensure Data sync ON, set name, click OK | n/a | FAIL | FAILED | `button[name="button-Save"]:visible` waitFor timed out 30s — no TableView ⇒ no SAVE button |
| 3 | Step 3: cancel auto-Share dialog | n/a | FAIL | FAILED | Cascade-failed (Step 2 didn't open the Save dialog) |
| 4 | Step 4: navigate to Browse > Dashboards, assert tile visible | n/a | FAIL | FAILED | Cascade-failed (no project saved). Tile name attr correctly uses verbatim `[name="div-UiSmoke...{stamp}"]` (PascalCase preserved per MCP observation) |
| 5 | Step 5-6: right-click tile → Share → fill recipient → OK | n/a | FAIL | FAILED | Cascade-failed |
| 6 | Step 7-8: double-click tile to reopen, verify demog table loaded | n/a | FAIL | FAILED | Cascade-failed |
| 7 | Step 9-11: right-click tile → Delete → confirm → assert tile gone | n/a | FAIL | FAILED | Cascade-failed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 3m 42s |
| **Total scenario run** | 3m 42s |

## Summary

Spec FAILED on dev.datagrok.ai (3m 42s). 7 soft-step(s) failed; root cause is Step 1 — the right-click → Open dispatch fires but does not materialize a TableView in the Playwright (`channel: 'chrome'`) runtime.

In a freshly-driven MCP session earlier today the SAME flow PASSED end-to-end (right-click → Open opens `/file/System.DemoFiles/demog.csv`, SAVE button visible, dialog opens with Data Sync ON, OK click triggers POST `/api/projects/{id}` → 200, auto-Share dialog appears, etc). The Playwright reproduction is the gap.

Hypothesis (pending new MCP session for live debugging):
- The synthetic `dispatchEvent('contextmenu')` on the `.d4-tree-view-node` may not deliver to Datagrok's Dart-side handler under Playwright the way it does under MCP (different CDP timing, different event bubbling chain).
- OR the `[name="div-Open"]` menu item click happens before the menu's own click handler is attached.
- OR the tree-node `:scope > .d4-tree-view-node` resolution picks an empty placeholder before lazy-load completes for the Demo subtree.

Next step (deferred to new MCP session per user direction): replay the right-click + Open click in MCP browser, then probe `grok.shell.tv?.dataFrame` immediately after, capture timing.

## Retrospective

### What worked well
- Anchored DOM walk through `.d4-tree-view-root > .d4-tree-view-group-host > children[3] (Files) > children[2] (Demo)` correctly locates the Demo subtree.
- `[name="div-Open"]` selector for the context-menu Open item is verified (per MCP observation block at top of spec).
- Tile slug rule fix: spec now uses `grok.shell.project.name` verbatim (no `.toLowerCase()`) to match `[name="div-UiSmoke<stamp>"]`.

### What did not work
- Synthetic context-menu dispatch + Open click does not materialize a TableView in Playwright `channel: 'chrome'` even though MCP-driven equivalent does.

### Suggestions for the platform
- The toolbar SAVE button's offsetWidth=0 behavior under JS-API openTableFromFile vs visible under right-click → Open (bug 2b) is the root differentiator that makes the UI Save flow brittle.

### Suggestions for the scenario
- Capture exact timing + event sequence in MCP for the right-click → Open path; mirror that in Playwright (likely needs to wait for `.d4-tree-view-node` to be a real node, not a placeholder, before dispatching).

## Raw error

```
Error: 7 step(s) failed:
  - Step 1: open demog.csv (UI — Browse-tree right-click → Open):
    expect(received).toBeGreaterThan(expected)
    Matcher error: received value must be a number or bigint
    Received has value: undefined

  - Step 2: open Save dialog, ensure Data sync ON, set name, click OK:
    locator.waitFor: Timeout 30000ms exceeded.
    waiting for locator('button[name="button-Save"]:visible').first() to be visible

  - Steps 3-9-11: cascade-failed (no project to operate on).
```
