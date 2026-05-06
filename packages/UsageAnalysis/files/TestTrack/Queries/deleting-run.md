# Queries — Deleting a Query — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse → Databases → Postgres → NorthwindTest, find `new_test_query` | 9s | PASS | PASSED | UI label `NorthwindTest` is the `friendlyName` of the `PostgresTest` connection (`a2d74603-…`). Double-clicked `Postgres` then `NorthwindTest`; `grok.shell.v.type` flipped to `queries`. The query was missing from prior runs — seeded via `conn.query('new_test_query', 'select * from orders')` + `grok.dapi.queries.save`. The card appeared in the gallery within ~600ms. |
| 2 | Right-click query → Delete → confirm DELETE | 6s | PASS | PASSED | First `dispatchEvent('contextmenu')` without `clientX/clientY` opened a popup container with no labels; re-dispatching with viewport coords surfaced the full menu (`Edit…, Clone, Run, Debug, Rename…, Clear cache, Delete, Share…, …`). Clicking `.d4-menu-item-label` "Delete" opened the confirm dialog (`Are you sure? Delete query "new_test_query"?`). The DELETE button has no `name=` and is not a direct child of `.d4-dialog`; it lives in `.d4-dialog-footer`. After click, `grok.dapi.queries.find(id)` returned null on the first poll (~300ms). |
| 3 | Refresh Browse — verify deletion | 2s | PASS | PASSED | The gallery card disappeared automatically after server-side delete; no explicit Refresh click needed. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 38s |
| grok-browser execution (scenario steps) | 17s |
| Execute via grok-browser (total) | 55s |
| Spec file generation | 30s |
| Spec script execution | 25s |
| **Total scenario run (with model)** | 1m 50s |

## Summary

Delete flow worked in both MCP and Playwright. The dev server's `NorthwindTest` is the
friendlyName of the `PostgresTest` Postgres connection; the canonical name is
`PostgresTest`. The new MCP run uncovered two refinements over the previous run:
the `contextmenu` event must include `clientX/clientY` for the menu's labels to render,
and the confirmation `DELETE` button is reliably reachable via `.d4-dialog-footer button`
rather than `.d4-dialog button`. **Total scenario run (with model): 1m 50s.**

## Retrospective

### What worked well
- `dispatchEvent('contextmenu')` with viewport coordinates (`clientX/clientY` from `getBoundingClientRect`) reliably opens the gallery context menu with all items rendered.
- Server-side deletion completes in <300ms; the gallery removes the card without an explicit Refresh.
- `.d4-dialog-footer button` filtered by `hasText: /^DELETE$/` is a stable Playwright selector for the confirmation button.

### What did not work
- A bare `dispatchEvent('contextmenu')` (no coords / no `buttons: 2`) opened a popup container but rendered zero `.d4-menu-item-label` children — the menu is built only when synthetic mouse coordinates are present. First attempt looked like a no-op.
- `.d4-dialog button` doesn't always select the visible footer buttons depending on dialog markup; scoping to `.d4-dialog-footer` is more robust.
- `grok.dapi.connections.filter('name = "NorthwindTest"').first()` returns `MSSQLTest` because the smart-search filter is substring-based and ranking-dependent — confusing when scripting setup. Always include `dataSource = "Postgres"` or look up by ID.

### Suggestions for the platform
- Make the delete confirmation dialog buttons carry `[name="button-DELETE"]` / `[name="button-CANCEL"]` attributes — eliminates text-matching fallback.
- The "right-click" path on a gallery card should not depend on the synthetic event's `clientX/clientY` to decide whether to render menu labels — either render unconditionally or document the expectation. This silently broke the first attempt.
- Consider exposing `grok.dapi.connections.byFriendlyName(name)` (or a strict `eq()` operator) so scripted setup doesn't return another connection that contains the substring.

### Suggestions for the scenario
- Step 2 bullet "Right-click the **connection**" should say "Right-click the **query card**". As written, readers would try deleting `NorthwindTest` itself.
- Step 3 says "Reftesh" — typo, should be "Refresh".
- Add a pre-condition: "The query `new_test_query` must exist (left behind by the edit scenario)." Without it, readers hit a missing-query error.
- Clarify that `NorthwindTest` is the friendlyName of the `PostgresTest` connection on dev — newcomers script-looking it up by canonical name will get a different connection.
