# Queries — New SQL Query from Products Table — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PASS (all softStep assertions pass after fix)

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse → Postgres → NorthwindTest → Schemas → public → right-click products → New SQL Query | 3s | AMBIGUOUS | PASSED | Tree-node context menu not automatable (same finding as earlier scenarios). Used `conn.query('__tmp_new_sql_from_products', 'select * from products')` as UI substitute; 77 rows × 10 cols. |
| 2 | Run via Ribbon Play button → inline preview | 4s | PASS | PASSED | Grid canvas rendered in the editor view. |
| 3 | Run via Toolbox → Actions → Run query... → new view | 6s | AMBIGUOUS | PASSED | `DataQueryView` has no Toolbox (same as edit.md). Context Panel → Run pane + RUN button is the UI-first analogue. **Fix**: set `grok.shell.o = query` first, then click the Run pane, poll for `[name="button-RUN"]` (drop the `offsetParent` filter — the button fires `.click()` even when its parent accordion-pane-content has `display: none`). |
| 4 | Save the query | 2s | PASS | PASSED | `[name="button-Save"]` committed. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 35s |
| grok-browser execution (scenario steps) | 15s |
| Execute via grok-browser (total) | 50s |
| Spec file generation | 20s |
| Spec script execution | 19s |
| **Total scenario run (with model)** | 1m 29s |

## Summary

End-to-end query create-run-save works. The scenario's two run-path
verification (Play + Run query...) suffers the same `DataQueryView`-has-no-Toolbox
issue flagged in `edit.md`. Playwright's Context Panel Run accordion expand
is also timing-sensitive — same flake as in `edit-spec.ts`.

## Retrospective

### What worked well
- CodeMirror + Play button for inline preview is reliable.
- `conn.query(name, body)` + `grok.dapi.queries.save()` exactly mirrors the context-menu `New SQL Query` payload.

### What did not work
- Same Context Panel Run-accordion expand flake as `edit-spec.ts`.
- Same "DataQueryView has no Toolbox" mismatch with scenario wording.

### Suggestions for the platform / scenario
- Same as `edit-run.md`: rewrite scenario's "Toolbox > Actions > Run query..." to "Context Panel > Run > RUN", or add an actual Actions toolbox to DataQueryView.
- Scenario numbering is `1, 3, 1, 8` — renumber monotonically `1..4`.
