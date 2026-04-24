# Queries — Deleting a Query — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse → Databases → Postgres → NorthwindTest and find `new_test_query` | 7s | PASS | PASSED | Navigated `/browse`, double-click Postgres and NorthwindTest; `new_test_query` visible in the 14-card queries gallery. |
| 2 | Right-click query → Delete → confirm DELETE | 8s | PASS | PASSED | `contextmenu` dispatched on the gallery card → menu showed 14 items (Edit, Clone, Run, Debug, Rename, Clear cache, `Delete`, Share, ...). Clicked `Delete`. A confirmation dialog appeared (`Are you sure? Delete query "new_test_query"?`) with `DELETE` / `CANCEL` buttons. Clicked `DELETE`. Verified via `grok.dapi.queries.find(id)` returning null. |
| 3 | Refresh Browse — verify deletion | 2s | PASS | PASSED | No element with text `new_test_query` visible in the gallery — the card disappeared without an explicit Refresh click. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 40s |
| grok-browser execution (scenario steps) | 17s |
| Execute via grok-browser (total) | 57s |
| Spec file generation | 35s |
| Spec script execution | 17s |
| **Total scenario run (with model)** | 1m 49s |

## Summary

The delete flow worked in both MCP and Playwright. Key finding: the scenario's
directive to right-click "the connection" is incorrect — the user must right-click
the **query card**, not the `NorthwindTest` connection. The confirmation dialog
uses text-matched DELETE/CANCEL buttons (no `name=` attribute) and the gallery
updates automatically after the server confirms.

## Retrospective

### What worked well
- `dispatchEvent('contextmenu')` on the gallery card opens the query context menu with all expected items.
- Delete → DELETE flow is one click each; server-side deletion is immediate (~300ms).
- Gallery auto-removes the deleted card without a manual Refresh click.

### What did not work
- **Scenario step 2 wording is wrong**: "Right-click the **connection**" — the connection deletion is a different flow. To delete a query you must right-click the **query card** inside the queries gallery.
- Initial attempt to enumerate confirmation-dialog buttons via generic `.d4-dialog` selector failed: the dialog uses a slightly different class tree; the a11y snapshot's `dialog` role was the reliable signal (button text `DELETE`).

### Suggestions for the platform
- Make the delete confirmation dialog buttons carry `[name="button-DELETE"]` / `[name="button-CANCEL"]` attributes like standard dialogs — it removes the text-matching fallback.
- Gallery refresh after delete is automatic — document that this is different from Browse-tree refresh (which does need an explicit click in some cases).

### Suggestions for the scenario
- Step 2 bullet "Right-click the connection" should say "Right-click the query card". As written, readers would try deleting `NorthwindTest` itself.
- Step 3 says "Reftesh" — typo, should be "Refresh".
- Add a pre-condition: "The query `new_test_query` must exist (left behind by the edit scenario)." Without the precondition readers will hit a missing-query error.
