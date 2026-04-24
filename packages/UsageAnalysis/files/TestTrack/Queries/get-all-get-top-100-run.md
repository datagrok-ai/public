# Queries — Get All / Get Top 100 — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

### Part 1 — PostgresDart
| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1.1 | Browse → PostgresDart → NorthwindTest → Schemas → public | 1s | FAIL | SKIPPED | No `NorthwindTest` under PostgresDart on dev. |
| 1.2–1.4 | Right-click orders → Get All / Get Top 100 | — | SKIP | SKIPPED | Blocked by 1.1. |

### Part 2 — Postgres
| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 2.1 | Browse → Postgres → Northwind → Schemas → public | 2s | AMBIGUOUS | — | Dev has `NorthwindTest` (not `Northwind`). |
| 2.2 | Right-click orders → Get All | 4s | AMBIGUOUS | PASSED (JS) | Table-node context-menu is not automation-friendly (same finding as columns-inspect). JS substitute: `select * from orders` → 830 rows × 14 cols. |
| 2.3 | Right-click orders → Get Top 100 | 3s | AMBIGUOUS | PASSED (JS) | JS substitute: `select * from orders limit 100` → 100 rows. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 35s |
| grok-browser execution (scenario steps) | 8s |
| Execute via grok-browser (total) | 43s |
| Spec file generation | 30s |
| Spec script execution | 9s |
| **Total scenario run (with model)** | 1m 22s |

## Summary

End-to-end query execution against `NorthwindTest.orders` works — 830 rows for
the full fetch and 100 rows for the limited fetch. The scenario's UI-first
path (right-click table → Get All / Get Top 100) could not be exercised
reliably under automation because the table tree node does not reliably
expand / accept `contextmenu` dispatch. Part 1 is blocked because
`NorthwindTest` does not exist under PostgresDart on dev.

## Retrospective

### What worked well
- Direct query execution via `conn.query(...).executeTable()` is a reliable JS substitute for `Get All` / `Get Top 100`.
- Row counts match expected semantics (100 vs all).

### What did not work
- **Right-click on a Schemas table node is not automatable**. Same finding as `columns-inspect.md`.
- **No `NorthwindTest` under PostgresDart** on dev.

### Suggestions for the platform
- Expose a DOM-addressable `Get All` / `Get Top 100` action on the table node (e.g. an icon button with `name="button-GetAll"`) so automation can trigger it without the context menu.
- Seed a Northwind database under PostgresDart on dev for parity with the Postgres provider.

### Suggestions for the scenario
- Replace "Northwind" with "NorthwindTest" and flag PostgresDart as optional if the dev seeding is incomplete.
- Add a clarifying sentence that these are context-menu actions on a table inside `Schemas → public`, not on the table view itself.
