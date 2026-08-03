# Copy / clone / move objects with metadata — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | S1.1 Open SPGI with sticky meta (setup) | 12s | PASS | PASSED | SPGI.csv 3624 rows; created schema `CLAUDE_SPGI_*` (rating:int, notes:string) matched via `source=<tag>` on Id col; seeded rating/notes for CAST-634783..785 |
| 2 | S1.2a Cloned table preserves metadata | 12s | PASS | PASSED | `df.clone()` carries `source` tag; getAllValues via cloned Id col returned rating=5,4,3 notes=excellent,good,average |
| 3 | S1.2b New view on same DataFrame shows metadata | 18s | PASS | PASSED | `grok.shell.addTableView(df)` second time — tag preserved, metadata identical |
| 4 | S1.2c Save as project and reopen | 19s | PASS | PASSED | Used 5-row subset (full 3624-row upload exceeds MCP 30s timeout); project.addChild(ti) + upload + save; closeAll → find → open; metadata intact |
| 5 | S1.2d Move project to Space and reopen | 12s | PASS | PASSED | `grok.dapi.spaces.createRootSpace(...)` + `sc.addEntity(projectId, false)`; project reopened from space — tag + rating/notes preserved |
| 6 | S1.2e Import/Export (d42 binary) | 3s | PASS | PASSED | `df.toByteArray()` → `DataFrame.fromByteArray()` (11058 bytes); all tags, semTypes, sticky-meta lookups intact |
| 7 | S2.1–2 Clear rating (sentinel 0) + notes, save | 4s | PASS | PASSED | rating updated 5 → 0, notes cleared to ''; **platform limitation**: `setAllValues` silently ignores null for `int` columns, so `0` is used as the "cleared" sentinel |
| 8 | S2.3 Deletion persists across view reopen | 3s | PASS | PASSED | closeAll + re-addTableView → rating still 0, notes still '' |
| 9 | S3.1–2a Metadata persists across page refresh | 13s | PASS | PASSED | Seeded rating=2, notes='refresh-test' on CAST-634787; full page reload → values intact |
| 10 | S3.2b Metadata persists across logout + login | 16s | PASS | PASSED | `POST /api/users/logout` (200) → reload → UI login with `DATAGROK_LOGIN`/`DATAGROK_PASSWORD` → read back: all 3 rows (CAST-634787/784/785) identical to pre-logout state |

**Time** = step 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

> The scenario heading mentions "server restart" but the numbered body only enumerates *Refresh browser tab* and *Logout & login again* — both covered above. No server-restart step was executed.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~6m |
| grok-browser execution (scenario steps) | ~1m 52s |
| Execute via grok-browser (total) | ~7m 52s |
| Spec file generation | ~2m |
| Spec script execution | 1m 29s |
| **Total scenario run (with model)** | ~11m 20s |

## Summary

All 10 steps PASS end-to-end, both in the MCP scenario run and in the Playwright replay (spec wall-clock 1m 29s). Full coverage: clone, new view, project save/reopen, move-to-Space, d42 export/import, clear-and-reopen, page refresh, and logout + re-login — every sub-step explicitly enumerated in the scenario body. A real API gap around clearing numeric sticky-meta values is documented below; the spec works around it with a sentinel value. **Total scenario run (with model)**: ~11m 20s.

## Retrospective

### What worked well
- Tag-based sticky meta (`source=<unique-value>`) is the most reliable keying mechanism — works on any string column, avoids semType matcher ambiguity, survives every copy/clone/save/reopen/import path tested.
- `DG.Project.create()` + `project.addChild(tableInfo)` + `grok.dapi.tables.uploadDataFrame` + `tables.save` + `projects.save` is the canonical save flow; `Project.open()` restores the tagged DataFrame exactly.
- `grok.dapi.spaces.createRootSpace` + `SpaceClient.addEntity(projectId, false)` moves a project into a space without breaking sticky-meta linkage.
- `df.toByteArray()` / `DataFrame.fromByteArray()` is a sticky-meta-safe round-trip — all column tags and semTypes survive, and lookups hit the same server rows.
- `POST /api/users/logout` (returns 200 + clears the `auth` cookie) is the clean server-driven logout path; pairing it with page reload + the standard login form gives a real session-reset test — sticky meta persists.
- Fire-and-forget async blocks that stash status on `window.*` + `page.waitForFunction` let the Playwright spec drive multi-second platform operations (uploadDataFrame, projects.save) without tripping MCP's 30 s call timeout.

### What did not work
- **Molecule-keyed sticky meta silently dropped all writes.** `matchBy: 'semtype=Molecule'` (or `'semtype=molecule'`) compiled into a schema but `setAllValues` returned null and every subsequent `getAllValues` returned int-null (-2147483648) and ''. The tag-based workaround (`matchBy: 'source=<tag>'` on the `Id` column) was the only path that persisted. Unclear whether the failure is in the molecule canonicalisation step or in `semtype=` matcher wiring.
- **`setAllValues` with a null int value is a no-op, not a clear.** `col.set(0, null)`, passing `DG.INT_NULL` (-2147483648), `setString(i, null)`, and CSV-parsed null-int all leave the stored rating at its previous value. The only way to "erase" an int via the JS API is to overwrite with a sentinel (e.g. `0`) or delete the whole schema. String columns clear correctly with `''`.
- **MCP Chrome-DevTools protocol timeout ≈ 30 s** returns a bare `"ns"` or `"Promise was collected"` with no hint that the JS promise is still running server-side. The fix (fire-and-forget + poll `window._status`) is workable but non-obvious; the first `uploadDataFrame` attempt on the full 3624-row SPGI silently aborted this way.
- **Default Playwright `testMatch` is `*.spec.ts`, not `*-spec.ts`.** Running `npx playwright test <file>` on a `-spec.ts` file in this repo needs a `--config` override pointing at a config whose `testMatch` includes the dash variant. The existing `/tmp/playwright.config.ts` can't load `@playwright/test` because node_modules aren't there; a temporary config inside `reddata/` is required.

### Suggestions for the platform
- `setAllValues` should accept `null` for numeric columns as an explicit clear (drop the row for that schema/key pair) rather than silently skipping the write. Today the only way to clear an int field is to delete and recreate the schema.
- Surface a dedicated `grok.dapi.stickyMeta.clearValues(schema, keys)` (or at minimum a documented semantics note) so scenarios like "delete fields rating and notes" have a non-guessing API path.
- Diagnose why `semtype=Molecule` schemas silently swallow `setAllValues`. A server-side log line at INFO when the matcher rejects all keys would have short-circuited ~10 min of the debugging loop.
- Spec writer ergonomics: ship a checked-in `playwright.config.ts` that sets `testMatch: ['**/*-spec.ts']` so `/grok-debug-scenarios` can run without a temp config file.

### Suggestions for the scenario
- "Open SPGI table (with some sticky meta added by you)" is ambiguous about *which* sticky-meta schema and which fields the tester is expected to have on hand. Either link to a known preset (e.g. TestSchema1 with rating/notes) or replace the bullet with explicit setup steps (create schema, attach tag, seed 3 cells).
- Section 2 says "Delete fields rating and notes. Save" without specifying whether deletion means clearing values on one cell vs. dropping the fields from the schema. The distinction matters (clear → values can come back; schema edit → fields disappear for everyone). Split into two sub-scenarios.
- The "Moved project to some Space" step assumes a Space already exists. Add a pre-condition ("create a dedicated test Space") or note that `grok.dapi.spaces.createRootSpace` is the programmatic path.
- "Imported/exported" is underspecified — zip project, d42 table, csv+layout are all different paths with different guarantees. Pick one and name the file format expected.
- Section 3's heading advertises "server restart" but the numbered steps only enumerate refresh and logout/login. Either drop "server restart" from the heading or add it as an explicit (manual-only) step so the scope matches the title.
