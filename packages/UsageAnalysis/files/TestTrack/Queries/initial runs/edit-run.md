# Queries — Editing a SQL Query — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Right-click test_query → Edit... (open editor) | 12s | PASS | PASSED | Navigated Browse → Databases → Postgres (dblclick) → NorthwindTest (dblclick) → queries gallery, then dispatched `contextmenu` with viewport `clientX/clientY` on `div[name="div-TestQuery"]`. The `Edit...` menu item (`.d4-menu-item-label`) opened `DataQueryView` with `input-Name="test_query"` and `CodeMirror="select * from products"`. Spec uses a direct URL (`/queries/Dbtests.PostgresTest?browse=db`) to skip the inconsistent tree-expand path. |
| 2 | Change name to new_test_query and click SAVE | 4s | PASS | PASSED | `click` → `Ctrl+A` → `type('new_test_query')` → `Tab` (Dart change listener requires real keyboard events). `[name="button-Save"]` (Title-case — no all-caps variant exists). Verified via `grok.dapi.queries.find`: `name=NewTestQuery`, `friendlyName=new_test_query`. |
| 3 | Change query body: select * from orders | 2s | PASS | PASSED | `document.querySelector('.CodeMirror').CodeMirror.setValue('select * from orders')`. UI typing into CodeMirror remains unreliable. |
| 4 | Run via Ribbon Play button → inline preview grid | 3s | PASS | PASSED | `[name="icon-play"]` click; 3 grid canvases rendered inside the editor view. |
| 5 | Run via Toolbox > Actions > Run query... → new view | 5s | PASS | PASSED | **DataQueryView DOES have a Toolbox `Actions > Run query...`** (rendered as `label.d4-link-action`) — this contradicts the previous edit-run.md (2026-04-24) which claimed the Toolbox was missing and had to fall back to Context Panel → Run accordion. Click opens new `TableView` named `new_test_query`; the editor's `.CodeMirror` element detaches from the DOM (`cmCount: 0`). Spec waits for the detach as a sync barrier — Dart caches the editor's pending body during the unmount, so the next remount restores `orders`. |
| 6 | Save the query | 4s | PASS | PASSED | Switched back via `grok.shell.v = qv`, then **polled `cm.getValue() === 'select * from orders'`** (up to 10s) so the SAVE click happens only after the editor has remounted with the cached body. Without this remount-wait, SAVE persists the prior server-side body. Verified `friendlyName=new_test_query`, `query=select * from orders`. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 35s |
| grok-browser execution (scenario steps) | 14s |
| Execute via grok-browser (total) | 1m 0s |
| Spec file generation | 40s |
| Spec script execution | 3m 0s |
| **Total scenario run (with model)** | 6m 30s |

Spec script execution covers five `npx playwright test` runs total: 1m 6s (initial — tree-expand failed, switched to direct-URL), 31s + 36s + 36s (three failed save-step patches), 21s (final pass after diagnosing the editor unmount/remount cycle).

## Summary

All six scenario steps pass in both the MCP run on dev (https://dev.datagrok.ai/) and
the Playwright spec (`edit-spec.ts`, 21s headed run). The previous 2026-04-24 run-md
note that `DataQueryView` has no Toolbox is wrong — `Toolbox > Actions > Run query...`
exists (`label.d4-link-action`) and was exercised directly. The interesting failure
mode discovered: clicking `Run query...` opens a new `TableView` and detaches the
editor's `.CodeMirror` from the DOM. If the spec switches back to the editor too
quickly, Dart re-mounts the editor by reading the last server-saved body
(`select * from products`), discarding the unsaved `orders` edit; the next SAVE
then persists `products`. The fix is two waits: (a) wait for `cmCount === 0` after
clicking `Run query...` so Dart has cached the editor state, and (b) after
`grok.shell.v = qv`, poll `cm.getValue() === 'select * from orders'` until the
remount has restored the cached body, then click SAVE. Final state:
`name=NewTestQuery`, `friendlyName=new_test_query`, `query=select * from orders`.
**Total scenario run (with model): ~6m 30s.**

## Retrospective

### What worked well
- **Editor unmount/remount as a sync barrier**: `Toolbox > Actions > Run query...` deterministically detaches `.CodeMirror` (`cmCount: 0`). Waiting for that detach (then for the post-switch remount with the cached body) is a reliable signal for "Dart finished caching the pending editor state."
- **`Toolbox > Actions > Run query...`** is present on `DataQueryView` (selector: `label.d4-link-action`). Spec exercises the scenario's literal wording.
- **Direct URL entry** to `/queries/Dbtests.PostgresTest?browse=db` lands on the gallery in one navigation, sidestepping the inconsistent Browse-tree expand semantics (`Databases` 1-click vs `Postgres` dblclick).
- **`contextmenu` with viewport `clientX/clientY` on the gallery card** reliably surfaces the full menu (including `Edit...`).
- **Keyboard-driven name edit** (`click → Ctrl+A → type → Tab`) committed `friendlyName=new_test_query`; programmatic `.value` would have been silently dropped.

### What did not work (and how we fixed it)
- **Naïve view-switch in step 6**: `grok.shell.v = qv` + 500 ms wait + click `button-Save` reliably saves the *prior* server-side body in headed Playwright. Root cause: when `Run query...` opens the new TableView, Dart unmounts the editor's CodeMirror; switching back too fast causes Dart to re-mount the editor by reloading the last persisted body, discarding the unsaved edit. **Fix**: poll until `cmCount === 0` (editor cached) after step 5, and until `cm.getValue() === 'select * from orders'` after the switch-back, before clicking SAVE.
- **Defensive `cm.setValue` and direct `qv.query.query =` writeback before SAVE did not help** — the SAVE handler appears to read from a Dart-side model that doesn't observe these JS-side mutations, so the unmount/remount-cache path is the only reliable way to keep the unsaved edit alive.
- **First-pass spec used live tree expansion** (`Databases` click → `Postgres` dblclick → `NorthwindTest` dblclick → wait for the gallery card) and timed out at the `Postgres` step. Replaced with the direct URL.

### Suggestions for the platform
- The editor's "cache pending body during unmount, restore on remount" semantics are subtle and timing-sensitive — automation-friendly hooks (e.g. an `onEditorPendingFlush` event, or persisting the pending body synchronously into `query.query`) would eliminate the need for the cmCount-based barrier.
- Add a `[name="button-RUN-Query"]` (or any stable `name=`) on the Toolbox `Run query...` link — currently it's only addressable by text on `label.d4-link-action`.
- Standardize `button-Save` casing — every other ribbon button (`button-OK`, `button-CANCEL`, `button-DELETE`) is all-caps; only `button-Save` is Title-case. Either flip to `button-SAVE` or document the exception.
- Browse-tree per-provider grouping (e.g. `Postgres`) requires `dblclick` while sibling levels (`Databases`) expand on single click — automate-friendly trees should be uniform.

### Suggestions for the scenario
- The previous run's note "DataQueryView has no Toolbox" is wrong. Step 5 wording (`Toolbox > Actions > Run query...`) IS achievable — keep the wording as-is.
- Step numbering still jumps `1, 2, 3, 1, 8` in the source file — should be `1..6`.
- Add an explicit pre-condition that `test_query` must already exist on `NorthwindTest` (currently only implicit via `adding.md`'s order=1).
- Note that `NorthwindTest` (Browse label) is the friendlyName of the underlying `PostgresTest` connection on dev; scripted lookups by canonical name will get a different connection.
- Note that the saved entity's `name` will be normalized to PascalCase (`NewTestQuery`) while `friendlyName` stays `new_test_query` — useful for testers verifying persistence.
- Add an explicit final verification: "Confirm the query body is `select * from orders` after re-opening (Browse → query → Edit)" so the step-6 SAVE is asserted, not just clicked.
