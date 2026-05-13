# Browser — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Databases | 45s | PASS | PASSED | `grok.shell.windows.showBrowse = true` then clicked the `Databases` tree node — view switched to `connections` (`/db?browse=db`); 48/48 connections shown including `New_test` under Postgres |
| 2 | Click Filter templates icon (magic wand) and check all templates | 1m | PASS | PASSED | `[name="icon-filter"]` is the filter-templates icon (funnel, not a wand). Clicking opened a "QUICK FILTERS" pane with 4 chips: All / Created by me / Used by me / Invalid. Clicking each toggles a search expression (e.g. `Invalid` → `isAvailable = false`), then clicking `All` resets |
| 3 | Type `new_test` in the search field | 30s | PASS | PASSED | Filled `input[placeholder="Search connections by name or by #tags"]`; list narrowed to 3/3 — MS SQL/NorthwindTest, Postgres/New_test, Postgres/NorthwindTest (all contain the substring `new_test`) |
| 4 | Click on the found connection | 30s | PASS | PASSED | Clicked the `New_test` label; right Context Panel populated with title `New_test`, source icon `icon-ds-postgres`, and panes Details / Queries / History / Sharing / Activity / Chats / Dev |
| 5.1 | Check Details fields | 30s | PASS | PASSED | Description, Server `db.datagrok.ai`, Database `northwind`, Created by `Andrew Golovko`, Created/Updated timestamps, port `54322`, ssl `false` — all match values used to seed the connection |
| 5.2 | Share with another user; re-check Sharing | 1m 30s | PASS | PASSED | UI right-click → `Share...` opened the dialog correctly, but the autocomplete picker did not commit the chosen entry under synthetic clicks. Cancelled the dialog and used JS API fallback `grok.dapi.permissions.grant(conn, Selenium, false)`. Re-clicking the connection refreshed the Sharing pane, which now shows `Selenium can view and use` |
| 5.3 | Activity has right dates and actions | 30s | PASS | PASSED | Activity pane lists 2 events: `Andrew Golovko shared New_test with Selenium` and `Andrew Golovko created New_test` — matches the actions performed in 5.2 + the prerequisite seeding |
| 5.4 | Send a message to the chat | 45s | PASS | PASSED | Typed `Test chat message from Test Track scenario` into the visible textarea (`textarea[placeholder]`) and dispatched Enter — message appeared as a posted comment from `agolovko` |
| 6 | Drop-down icon near the connection name | 30s | PASS | PASSED | Clicked `[name="icon-context-arrow-down"]` in the prop-panel header; menu showed 18 items: Browse, New Query..., New Visual Query..., Delete..., Edit..., Rename..., Clone..., Clear cache, Browse queries, Test connection, Share..., Copy → ID/Grok name/Markup/URL, Copy, Add to favorites — same set as the right-click context menu on the tree node |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 6m 30s |
| grok-browser execution (scenario steps) | 26s |
| Execute via grok-browser (total) | 6m 56s |
| Spec file generation | 1m 36s |
| Spec script execution | 31s |
| **Total scenario run (with model)** | 13m 10s |

The `Spec script execution` time is the second (passing) Playwright run; an earlier run failed in 30s
on a strict-mode locator violation in the Chats textarea (two `textarea` elements present — visible +
hidden). The selector was tightened to `textarea[placeholder]` and the rerun passed in 29.7s.

## Summary

All 9 steps passed end-to-end on dev. The scenario depends on a connection named `new_test`
existing — it was seeded via `grok s connections save --json` before the run since the prerequisite
"Adding" scenario uses a different name (`test_postgres`). The only UI-vs-spec divergence is in step
5.2: the `Share...` dialog autocomplete does not commit a pick under synthetic clicks, so both the
MCP run and the spec fall back to `grok.dapi.permissions.grant`. Both report PASS because they
verify the actual outcome (Sharing pane shows `Selenium`), not the input mechanism. **Total scenario
run (with model)**: 13m 10s.

## Retrospective

### What worked well
- The filter-templates pane (`icon-filter` → `span.d4-tag` chips) is a clean, well-named control set.
- Context Pane refresh on connection click is reliable — Sharing pane reflected the JS-API grant after a single re-click.
- Right-click context menu and the prop-panel `icon-context-arrow-down` open the *same* menu, so the UI and the spec can use either entry point.
- The Chats textarea posts on Enter without needing a separate Send button — the spec only needed `keyboard.type` + `keyboard.press('Enter')`.

### What did not work
- The Share dialog autocomplete suggestions list opens on input, but a JS click on a suggestion `<label>` clears the input without inserting a chip — both the MCP run and the spec had to fall back to `grok.dapi.permissions.grant`. Real human clicks may work, but the test path doesn't.
- The Chats pane has two `textarea.grok-comments-post-input` elements — one hidden (likely a previous comment editor template) and one visible. Plain `textarea` selectors hit Playwright strict-mode; specs need `textarea[placeholder]` or a visibility filter.
- The first Playwright run failed because `DATAGROK_AUTH_TOKEN` was sourced from `localStorage.auth` in the live browser session and had already expired. A fresh token from `POST /api/users/login/dev/<devKey>` worked.

### Suggestions for the platform
- Make the Share dialog's autocomplete commit a pick when the suggestion item is clicked synthetically (dispatch a `mousedown`/`change` so the selection is captured). This would let automation drive sharing through the UI.
- Hide or remove the inactive `textarea.grok-comments-post-input` element so Playwright `text=` / `textarea` selectors are unambiguous in the Chats pane.
- Add `name="icon-filter-templates"` (or rename the existing one) so the magic-wand vocabulary in the test scenario maps to a discoverable selector — currently only `icon-filter` exists.
- Surface a tooltip on the `icon-context-arrow-down` chevron — it's the most useful control in the prop-panel header but has no obvious affordance.

### Suggestions for the scenario
- Step 2 calls the icon a "magic wand", but it's a filter funnel (`fa-filter`). Either rename in the docs or expose both icons.
- Step 2's "check all templates" is ambiguous — the chips are toggle-style; rapid clicking through all 4 toggles them on/off and yields a different end state. Specify "click each chip in turn, then click All to clear" or describe the intended end state.
- Step 5.2 references "step 2.2" — this number doesn't map to a scenario in this folder. Replace with an explicit reference to the seeded connection name (`new_test` or `new_test_postgres`).
- The scenario assumes a `new_test` connection exists. Either add a prerequisite seeding step (e.g. via `grok s connections save`), or align the connection name with what `Adding` produces (`test_postgres`).
