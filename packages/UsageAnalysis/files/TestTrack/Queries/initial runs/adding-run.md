# Queries — Adding a New SQL Query — Run Results

**Date**: 2026-05-04
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse → Databases → Postgres | 8s | PASS | PASSED | `Databases` expanded on single click; `Postgres` required `dblclick` event to expand. `NorthwindTest` visible after ~3.6s tree wait. |
| 2 | Right-click NorthwindTest → New Query... | 4s | PASS | PASSED | Dispatched `contextmenu` on `.d4-tree-view-node`; clicked `New Query...` menu item by text. URL updated to `/query/{uuid}`, page title `New query`, CodeMirror + `input[name="input-Name"]` mounted. |
| 3 | Enter `test_query` into Name field | 3s | PASS | PASSED | Located by `input[name="input-Name"]` directly (simpler than the prior `offsetParent + placeholder` heuristic). Click + Ctrl+A + `pressSequentially` committed value as `test_query`. |
| 4 | Enter SQL `select * from products` | 2s | PASS | PASSED | JS API fallback — `document.querySelector('.CodeMirror').CodeMirror.setValue(...)`. UI typing into CodeMirror is unreliable. |
| 5 | Run via Ribbon Play button (inline preview) | 4s | PASS | PASSED | Clicked `[name="icon-play"]`; 3 grid canvases rendered inside the query view immediately. |
| 6 | Run via Toolbox → Actions → Run query... (new view) | 5s | PASS | PASSED | Clicked the `Run query…` toolbox label; new `TableView` named `test_query` opened in ~530 ms (`viewsBefore=2 → viewsAfter=3`). |
| 7 | Save the query | 8s | PASS | PASSED | Switched focus to `DataQueryView` (the inline run had grabbed focus), clicked `[name="button-Save"]`. Saved entity: `name=TestQuery`, `friendlyName=test_query`, `connection.name=PostgresTest`, `query=select * from products`. |

**Time** = step 2b wall-clock incl. model thinking. **Result** = MCP outcome. **Playwright** = spec step outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 25s |
| grok-browser execution (scenario steps) | 9s |
| Execute via grok-browser (total) | 34s |
| Spec file generation | 50s |
| Spec script execution | 21s |
| **Total scenario run (with model)** | 1m 45s |

## Summary

All seven scenario steps passed on dev (https://dev.datagrok.ai/) in both the MCP run
and the regenerated Playwright spec (`adding-spec.ts`, 17.6s headed run). The query was
created, executed via both run paths (inline preview + new-view), and persisted.
**Total scenario run (with model): ~1m 45s.** Persisted entity confirmed: `name=TestQuery`
(server normalizes to PascalCase), `friendlyName=test_query`, `connection=PostgresTest`,
`query=select * from products`. No regressions versus the prior 2026-04-24 run.

## Retrospective

### What worked well
- **`input[name="input-Name"]` directly addressable** — no need for the previous run's `offsetParent + placeholder` heuristic. Simplifies the spec.
- **`dispatchEvent('contextmenu')` on `.d4-tree-view-node`** opened the Browse context menu reliably (no native right-click needed).
- **CodeMirror `.CodeMirror.setValue(...)`** set the SQL body in a single call; UI keystrokes into CodeMirror remain unreliable.
- **Both run paths produced distinct surfaces** — Play button rendered an inline grid in the same view (3 canvases), Toolbox `Run query…` opened a separate `TableView` named `test_query`.
- **Playwright MCP fallback worked** — chrome-devtools MCP was stuck (`No page selected` on every call including `list_pages`); switching to Playwright MCP unblocked the run with no scenario changes.

### What did not work
- **chrome-devtools MCP server returned `No page selected` for every call** including `list_pages` and `new_page`, even after creating a tab via the raw DevTools `/json/new` HTTP endpoint. The MCP server appears to be holding a stale/null page handle. Recovered by using Playwright MCP instead.
- **Save button selector `[name="button-SAVE"]` does not exist** — actual selector is `[name="button-Save"]` (Title-case). The previous spec hedged with `[name="button-SAVE"], [name="button-Save"]`; the all-caps form has never matched on the query editor ribbon. Cleaned up to a single `[name="button-Save"]`.
- **`Postgres` tree node still requires double-click to expand** while `Databases` expands on single click — same inconsistency observed on 2026-04-24.

### Suggestions for the platform
- Make tree-node expand behavior uniform — either all levels expand on single click, or none do. The current `Databases` (1-click) vs. `Postgres` (dblclick) split is unique to the per-provider grouping layer and trips up automation each run.
- Standardize ribbon-button `name=` casing — the project uses `button-OK`, `button-CANCEL`, `button-DELETE` (all-caps), but `button-Save` (Title-case) breaks the pattern. Either normalize to all-caps or document the exception in [dialogs-menus.md](.claude/skills/grok-browser/references/dialogs-menus.md).
- Surface a visible indicator when a Name field's value will be normalized server-side (e.g. show `→ TestQuery` next to the input). Today the input keeps showing `test_query` while the saved entity is named `TestQuery`.
- Consider showing the underlying connection name in the Browse tree (or in the query editor header) when a display alias is in effect — `NorthwindTest` (Browse label) → `PostgresTest` (connection.name) is silent.

### Suggestions for the scenario
- The scenario numbering is broken: source lines are numbered `1, 3, 4, 5, 1, 8` — should be a monotonic `1..7`.
- Step 3 should specify the expected saved name (and whether normalization is intentional), so the tester knows whether `TestQuery` is a pass or a fail.
- Add an explicit final verification step: "Confirm the query appears under NorthwindTest in the Browse tree (refresh if needed)." The scenario currently stops at "Save", which does not assert persistence.
- Pre-condition: add "Ensure no query named `test_query` already exists on this connection" so re-runs don't collide with stale data.
- Step 1 lists `Browse > Databases > Postgres` but does not explicitly mention double-click for `Postgres`; reflect the inconsistent expand behavior or remove it as a platform fix.
