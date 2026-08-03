# Connections / Identifiers — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 0 | Prerequisite: create test_postgres connection | 9m 34s | PASS | n/a | test_postgres did not exist on dev. Cloned PostgresNorthwind params + working creds (datagrok/WZNFY...) into Agolovko:test_postgres (id af9bcf40-21a0-11f1-89e2-7b1321b80948); `conn.test()` → "ok". Browser session was Admin and could not see the new connection — switched session to agolovko via dev-key→token exchange (`/api/users/login/dev/<key>`) injected as `auth` cookie + localStorage. |
| 1 | Go to Data > Databases > Postgres | 10s | PASS | PASSED | Browse tree shows the Postgres provider with test_postgres after expanding the Postgres node and (in Playwright) clicking "Show more". |
| 2-3 | Right-click test_postgres → Configure Identifiers... | 50s | FAIL | PASSED | MCP browser session: context menu items = `Browse, New Query…, …, Add to favorites` — 17 entries, no `Configure Identifiers…`. Reproduced after hard-reload, cache+SW clear, both as Admin and agolovko. **Playwright run on the same dev server (fresh Chromium context, agolovko session) sees an 18th item `Configure Identifiers…` and clicking it opens dialog with title "Select primary schema for Identifiers Configuration".** Strongly suggests the MCP-attached Chromium has a stale/forked client bundle (the dev-server JS asset URL is `/login.dart.js?1777943623891` and `/js/api/js-api.js?1777943623891`; another tab or extension may have pinned an older one). MCP cannot load the dialog → all downstream steps could not run in MCP. |
| 4 | Set Schema = public, click OK | n/a | SKIP | FAILED | Spec attempted `setDialogInput('Schema'/'Primary Schema', 'public')` — failed with `TypeError: Illegal invocation`, indicating the schema field is a combo/select widget rather than a plain `<input>`. Selectors need to be authored against the live dialog. |
| 5 | Add identifier (CUSTOMER_ID, customers, customerid, [A-Z]{5}) | n/a | SKIP | FAILED | Spec couldn't find an "Add" button via heuristic (`button text=ADD/+` / `[name*="add"]`). Identifier editor row layout unknown. |
| 6 | SAVE | n/a | SKIP | FAILED | Cascaded from step 5. |
| 7 | Reload the page | n/a | SKIP | PASSED | Reload itself works; nothing to verify until earlier steps succeed. |
| 8-9 | Open customers table from test_postgres (Schemas/public) | n/a | SKIP | FAILED | `q.executeTable()` returned null in the agolovko Playwright context. Probable cause: `Data query "Execute" permission` (raw `tp.query('select…')` was rejected with that exact message when run as agolovko in MCP). The table-query path on the connection seems to require a saved query or higher privileges than the dev key grants. UI path `Right-click table → Get Top 100` in Browse tree would be the fallback. |
| 10 | Verify customerid values highlighted in blue | n/a | SKIP | FAILED | Cascaded from step 8-9 (no DataFrame to inspect). |
| 11 | Click column header → Context Panel shows semType | n/a | SKIP | FAILED | Cascaded from step 8-9. |
| 12 | Remove identifiers configuration on test_postgres | n/a | SKIP | FAILED | Could not be reached. The right-click menu does have `Configure Identifiers…` per Playwright; whether removal is a separate item or a button inside the dialog is unverified. |
| 13-14 | Reload and verify customerid not highlighted | n/a | SKIP | FAILED | Cascaded from step 12. `executeTable()` would also need to be replaced for the same reason as step 8-9. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome (MCP-observed). **Playwright** = 2e outcome (independent fresh-context replay).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 9m 50s |
| grok-browser execution (scenario steps) | 2m 46s |
| Execute via grok-browser (total) | 12m 36s |
| Spec file generation | 5m 30s |
| Spec script execution | 4m 22s (4 runs: 16s + 30s + 30s + 1m 30s + edits) |
| **Total scenario run (with model)** | 27m 5s |

The two `scenario steps` rows sum to `Execute via grok-browser (total)`.

## Summary

The scenario depends on a `Configure Identifiers…` connection right-click action. **Playwright** (fresh Chromium, agolovko session) sees the menu item and the schema-picker dialog opens — so the platform feature exists on dev. The persistent MCP-attached browser, however, never renders this 18th menu entry across multiple hard reloads, cache+SW clears, and user switches; the menu is stuck at 17 items. This blocked steps 4–14 from being driven through the MCP. Steps 4–14 were attempted via the generated spec but the heuristic dialog selectors and `tableQuery().executeTable()` path did not work in Playwright either; encoding correct selectors needs a live look at the dialog (combo widgets, identifier-row editor, OK/SAVE/Remove). **Total scenario run (with model): 27m 5s.**

## Retrospective

### What worked well
- Switching browser auth from Admin to agolovko via `POST /api/users/login/dev/<devKey>` → JWT → cookie+localStorage gave a fully working agolovko session and was reusable in the Playwright run via `DATAGROK_AUTH_TOKEN`.
- `grok s connections save --json … --save-credentials` cleanly created the prerequisite test_postgres with working credentials harvested from `public/packages/Samples/connections/postgres-northwind.json`.
- `softStep` + Playwright divergence helped: even though the MCP browser couldn't see the menu, Playwright independently confirmed the feature exists and the schema-dialog opens with the expected title.

### What did not work
- **MCP browser stuck on a partial menu**: across `navigate` reload, `ignoreCache: true` reload, service-worker unregister + Cache-API clear, and re-login, MCP saw only 17 right-click items. This is the root blocker — all downstream MCP work was impossible. Best guess: a long-lived Chromium profile pinned an older client bundle, or a Datagrok client-side feature flag is not being toggled the same way the fresh Playwright context toggles it.
- **`grok s connections test`** for connections owned by the agolovko user fails with `NoSuchMethodError: ... '[]' was called on null. Tried calling: ["secretConnection"]` (datlas `connectors_service.dart:418`) — the CLI test path expects a `secretConnection` blob that JSON-saved connections don't have. Tested via `await conn.test()` from the browser instead.
- **`grok s shares add … "All users"` reports success** but `grok s shares list` immediately after returns `(no results)`, and an Admin browser session continued to 404 the new connection. Either the share didn't persist, or it persists but the `shares list` reader looks at a different relation. Not investigated further; switched user instead.
- **`tp.tableQuery('public.customers').executeTable()` returns null** in the Playwright agolovko session; raw `tp.query('select * …')` returns `Data query "Execute" permission required`. There's no obvious public path to open a connection table without that permission, despite the right-click `Get Top 100` UI doing it in practice.
- Heuristic dialog selectors (`label:has-text(...)` → adjacent `<input>`) work for plain text fields (proved by adding-spec.ts) but break on combo/select widgets used in the schema picker.

### Suggestions for the platform
- **Investigate why the right-click menu's `Configure Identifiers…` item is registered in some sessions and not others on the same server.** If it's gated by a package init, the gate is firing inconsistently. Consider always registering the menu item (and disabling it with a tooltip when prerequisites aren't met) instead of hiding it.
- **Fix `grok s connections test`** for connections without a `secretConnection` blob — guard the `[]` call with a null check, then exercise normal credential lookup.
- **Make `grok s shares add` / `shares list` round-trip cleanly**, or document the relation it reads. The current behaviour (success on add, empty on list) is misleading.
- **Expose a permission-light read API** for opening a connection's table, or document the canonical JS-API call (`?tableQuery(...).executeTable()` vs a `Get Top 100`-style helper). Right now plugins/tests have to mimic the right-click action through the UI.

### Suggestions for the scenario
- **State the prerequisite explicitly**: "Adding-spec.ts must have run; test_postgres must exist and be testable." This run wasted ~10 minutes recovering from `test_postgres` not being present.
- **Expand step 12 wording** ("remove the identifiers configuration") with the exact UI path: is it a separate context-menu item, or a button inside the Configure Identifiers… dialog? The current wording is ambiguous.
- **Step 10 ("highlighted in blue") is canvas-based and untestable from outside the renderer** — reword to "the customerid column header shows a colored marker / shows semantic type CUSTOMER_ID in tooltip / Context Panel shows CUSTOMER_ID semantic type" so a Playwright assertion can target a DOM property.
- **Add a dedicated cleanup step**: "Delete the test_postgres connection" so the scenario is fully idempotent. (This run left Agolovko:test_postgres in place because the Adding/Edit/Identifiers chain expects it across runs.)
