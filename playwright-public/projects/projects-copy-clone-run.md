# Projects / Copy Clone: 3 save modes + GROK-19750 invariant — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Setup + 4b/c/d: open demog, add viewers, save 4 variants in one session | n/a | PASS | PASSED | tv-settle poll budget extended 12s → 30s (full-suite-run flake fix) |
| 2 | Step 4b verification: reopen `<name>-link`, table re-materializes (REOPEN #1) | n/a | PASS | PASSED | rowCount > 0 — `df.clone()` per save delivers independent reopen |
| 3 | Step 4b GROK-19750 invariant: reopen original, table re-materializes (REOPEN #2) | n/a | PASS | PASSED | rowCount > 0 — invariant holds |
| 4 | Step 5: re-share each variant via JS API | n/a | PASS | PASSED | Recipient = `target.group` (re-fetched via `dapi.users.find(id)` for materialized .group) |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 02s |
| **Total scenario run** | 1m 02s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 02s.

Three fixes converged this session:
1. `df.clone()` per save → each of the 4 variants (original / link / clone / pvc) gets an independent `tableInfo` (avoids the shared-tableInfo race that produced rowCount=0 on reopen of earlier variants).
2. `permissions.grant` recipient is `target.group` (Group, not User) — re-fetched via `dapi.users.find(id)` to ensure `.group` is materialized; avoids the `permissions_user_group_id_fkey` FK violation.
3. tv-settle poll budget extended from 12s → 30s — the full-suite run hit a 12s timeout under load while solo runs did not.

GROK-19750 invariant holds: after sibling Save Copies are saved, the original project still reopens with its data intact.

## Retrospective

### What worked well
- `df.clone()` per save isolates table state across variants.
- 30s tv-settle is sufficient under full-suite load.
- target.group recipient pattern avoids the FK violation cleanly.

### What did not work
- Nothing notable in the final run.

### Suggestions for the platform
- After `addTableView(...)` settles, `grok.shell.tv` should reliably expose `addViewer` within a tighter window than 30s.

### Suggestions for the scenario
- None from this run.
