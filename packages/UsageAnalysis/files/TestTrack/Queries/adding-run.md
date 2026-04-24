# Queries — Adding a New SQL Query — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse → Databases → Postgres | 12s | PASS | PASSED | Single-click expanded `Databases`; `Postgres` required `dblclick` to expand. `NorthwindTest` visible afterwards. |
| 2 | Right-click NorthwindTest → New Query... | 8s | PASS | PASSED | Dispatched `contextmenu` on `.d4-tree-view-node`; clicked menu item by text. Query editor (`.CodeMirror`) appeared, URL updated to `/query/{uuid}`. |
| 3 | Enter `test_query` into Name field | 6s | PASS | PASSED | Focused visible Name input, Ctrl+A + typed. Platform later normalized saved query name to `TestQuery` (PascalCase). |
| 4 | Enter SQL `select * from products` | 3s | PASS | PASSED | JS API fallback — `document.querySelector('.CodeMirror').CodeMirror.setValue(...)`. UI typing into CodeMirror is unreliable. |
| 5 | Run via Ribbon Play button (inline preview) | 5s | PASS | PASSED | Clicked `[name="icon-play"]`; grid canvas rendered inside the query view in ~600 ms. |
| 6 | Run via Toolbox → Actions → Run query... (new view) | 7s | PASS | PASSED | Clicked the `Run query...` label in Actions; new `TableView` named `test_query` appeared within ~850 ms. |
| 7 | Save the query | 17s | PASS | PASSED | Clicked `[name="button-SAVE"]`; record persisted. Stored connection is `PostgresTest` (same underlying JDBC connection the `NorthwindTest` Browse node points to), stored name normalized to `TestQuery`. |

**Time** = step 2b wall-clock incl. model thinking. **Result** = MCP outcome. **Playwright** = spec step outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 42s |
| grok-browser execution (scenario steps) | 16s |
| Execute via grok-browser (total) | 58s |
| Spec file generation | 1m 10s |
| Spec script execution | 27s |
| **Total scenario run (with model)** | 2m 35s |

## Summary

All seven scenario steps passed on the dev server in both the MCP run and the generated
Playwright spec (`adding-spec.ts`, 23.5s headed run). The query was created, executed via both
run paths (inline preview + new-view), and persisted. **Total scenario run (with model): ~2m 35s.**
Only one surprise: the stored entity was normalized to `TestQuery` under connection
`PostgresTest` — the `NorthwindTest` Browse node is a display alias over the `PostgresTest`
Postgres connection on dev. The scenario intent (create + run + save) succeeded.

## Retrospective

### What worked well
- `dispatchEvent('contextmenu')` on `.d4-tree-view-node` opened the Browse context menu reliably.
- CodeMirror `.CodeMirror.setValue(...)` fallback set the SQL body in a single call.
- Both run paths (Play button + Toolbox Actions → Run query...) produced distinct result surfaces (inline grid vs new view) as documented in `references/queries.md`.

### What did not work
- **`Postgres` tree node required a double-click to expand** — a single click selected but did not expand the subtree. `Databases` expanded on a single click. Inconsistent expand semantics between tree levels.
- **Save verification needed a fallback filter** — first attempt queried by `connection.name = "NorthwindTest"` but the saved query's connection is `PostgresTest`. The `NorthwindTest` Browse node label and the underlying connection entity name differ.
- **Query name normalization is silent** — `test_query` typed into the Name field is persisted as `TestQuery`. The input field value stays `test_query` visually, so this is not obvious at save time.

### Suggestions for the platform
- Make tree-node expand behavior uniform — either all levels expand on single click, or none do. Today the rule is inconsistent (Databases yes, Postgres no).
- Surface a visible indicator when a Name field's value will be normalized server-side (e.g. show "→ TestQuery" next to the input). Today users get a `test_query`-shaped input but a `TestQuery`-shaped entity.
- Consider displaying the underlying connection name in the Browse tree (or in the query editor header) when a display alias is in effect.

### Suggestions for the scenario
- The scenario numbering is broken: lines are numbered `1, 3, 4, 5, 1, 8` in the source file. Should be a monotonic `1..7`.
- Step 3 should specify the expected saved name (and whether normalization is intentional), so the tester knows whether `TestQuery` is a pass or a fail.
- Add an explicit final verification step: "Confirm the query appears under NorthwindTest in the Browse tree (refresh if needed)." — currently the scenario stops at "Save", which does not assert persistence.
- Pre-condition: add "Ensure no query named `test_query` already exists on this connection" so re-runs don't collide with stale data.
