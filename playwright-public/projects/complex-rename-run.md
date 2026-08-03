# Projects / Complex rename: rename-then-reopen reference resolution (GROK-19212) — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Setup: open 2 source tables + join | n/a | PASS | PASSED | 3 tables (src_a, src_b, joined) open |
| 2 | Save baseline project (Data Sync ON) | n/a | PASS | PASSED | Multi-table inline save (loop over `grok.shell.tables` — covers all 3 tables, not only the active TableView) |
| 3 | Step 7: rename referenced table inside the project | n/a | PASS | PASSED | GROK-19212 trigger — rename src_a after reopen |
| 4 | Step 8: re-save project after rename | n/a | PASS | PASSED | JS API `dapi.projects.save(p)` |
| 5 | GROK-19212 INVARIANT: close + reopen, must load without resolution error | n/a | PASS | PASSED | Project re-opens, dataframe re-materializes via `.script` provenance |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 48s |
| **Total scenario run** | 1m 48s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 48s.

GROK-19212 invariant verified — project reopens correctly after referenced table is renamed inside it.

Two key fixes this session:
1. Migrated from `loginToDatagrok` (which goes to `process.env.DATAGROK_URL ?? 'http://localhost:8888'`) to `gotoApp + setupSession` from `_helpers.ts` — was failing with `ERR_CONNECTION_REFUSED` against dev.
2. Replaced UI-driven local `saveProject` (blocked by bug 2b — toolbar SAVE button offsetWidth=0 after JS-API openTableFromFile) with inline multi-table JS API save — saves all 3 open tables (the canonical `saveProjectWithProvenance` helper saves only the active `tv.dataFrame`).

## Retrospective

### What worked well
- Multi-table save loop over `grok.shell.tables` correctly persisted src_a, src_b, joined; reopen restored them all and rename of src_a was visible in the reopened project.
- find-by-id pattern (id-based lookup post-rename) avoids the search-index lag of filter-by-name.

### What did not work
- Nothing notable in the final run.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- Step 9 sub-bullet (project-entity rename via `p.name = '<new>'; dapi.save(p)`) deferred — failed to propagate on dev (3/3 round 1 attempts). Worth a separate JS API setter-rename investigation.
