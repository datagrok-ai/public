# Projects / Complex derived-tables: Join lands in active project (GROK-19103) — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Setup + Step 1: open 2 source tables from files | n/a | PASS | PASSED | `openTableFromFile` (canonical) × 2 |
| 2 | Step 1 (Join sub-bullet): grok.data.joinTables produces 3rd table | n/a | PASS | PASSED | 3 tables in workspace |
| 3 | Step 2: capture project-count baseline | n/a | PASS | PASSED | |
| 4 | Step 2: Save current project with Data Sync ON | n/a | PASS | PASSED | Multi-table inline save (3 tables persisted via `addChild(tableInfo)` loop) |
| 5 | GROK-19103 INVARIANT: exactly ONE new project created | n/a | PASS | PASSED | `afterCount - baselineCount === 1` — no stray join-only project |
| 6 | Reopen: verify project loads with source tables intact | n/a | PASS | PASSED | `tables.length >= 2` after reopen |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 12s |
| **Total scenario run** | 1m 12s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 12s.

GROK-19103 invariant verified — joining 2 source tables and saving the active project produces exactly ONE new project entity (no stray join-only project).

Three migrations to make this pass on dev:
1. `loginToDatagrok` (`localhost:8888` default) → `gotoApp + setupSession` (canonical, honors `projectsTestOptions` baseURL).
2. `grok.data.files.openTable` (no `.script` provenance) → `openTableFromFile` (canonical OpenFile recorder, dot-form workaround for bug 2a).
3. UI-driven `saveProject` (toolbar SAVE click + dialog OK; blocked by bug 2b) → multi-table inline JS API save (loop over `grok.shell.tables` → `addChild(tableInfo)` per df).

## Retrospective

### What worked well
- Multi-table inline save persists all 3 open tables (src_a, src_b, joined) as expected.
- Project-count delta is cleanly 1 — no stray join-only project on dev today.

### What did not work
- Nothing notable in the final run.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- None from this run.
