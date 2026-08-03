# Projects / Lifecycle Files: open → save with provenance → reopen → share → rename — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 1m 00s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 00s |
| **Total scenario run** | 1m 00s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 00s.

Lifecycle Files: open via `openTableFromFile` → save with provenance via `saveProjectWithProvenance` → reopen → share → rename.

Fix this session: relaxed `.script` provenance regex from colon-only `OpenFile\("System:DemoFiles\/demog\.csv"\)` to `OpenFile\("System[:.]DemoFiles\/demog\.csv"\)` — accepts both colon-form and dot-form. The Phase B `openTableFromFile` helper normalizes `:` → `.` before passing fullPath to OpenFile (workaround for bug 2a, the colon-form Save POST silent-failure on dev).

Note: share step on this spec emits `Share skipped: ... permissions_user_group_id_fkey` (FK violation when granting to a User without materialized .group); soft-handled (warn, not assert) per spec design.

## Retrospective

### What worked well
- `openTableFromFile` produces `.script: 'Demog = OpenFile("System.DemoFiles/demog.csv") //...'` — reopen re-executes OpenFile and re-materializes the table.
- `saveProjectWithProvenance` mirrors the canonical `UITests/src/gui/gui-utils.ts:uploadProject` pattern.

### What did not work
- Nothing notable.

### Suggestions for the platform
- Per `PLATFORM-BUGS-FOUND-PHASE-B.md` #2a: platform should normalize the path-string in the Save dialog's `.script` serialization, or accept both forms cleanly so specs don't have to dance around the helper-side workaround.

### Suggestions for the scenario
- None from this run.
