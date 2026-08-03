# Projects / Lifecycle Spaces: GROK-18345 share + datasync — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: SKIP (env-blocked)

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Spaces lifecycle env-skip — addEntity probe | n/a | SKIP | SKIPPED | `addEntity` failed: `Operation caused an exception (pgpool1 ERROR 22P02 invalid input syntax for type uuid: "demog.csv")` |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | n/a (env-skipped) |
| **Total scenario run** | n/a |

## Summary

Spec env-skipped on dev.datagrok.ai due to bug 1 (Spaces `addEntity` UUID parse failure — see `PLATFORM-BUGS-FOUND-PHASE-B.md` #1).

The spec contains an up-front probe that calls `client.addEntity(saved, true)` with a saved FileInfo. The server-side payload encoding tries to parse `'demog.csv'` (the file name) as a UUID and fails with Postgres FK type-cast error. The probe catches this and calls `test.skip()` with a clear blocker message — the test will auto-unblock when the platform binding is fixed.

GROK-18345 invariant ("Recipient cannot open shared project that uses a Spaces dataset saved with data sync") cannot be exercised end-to-end at the JS API layer until the addEntity binding is fixed.

## Retrospective

### What worked well
- Probe-based env-skip pattern is robust — the test fails cleanly with a single skip message rather than producing 6 cascade-failed steps.

### What did not work
- Spaces lifecycle remains uncovered until platform-side fix.

### Suggestions for the platform
- Fix `dapi.spaces.id(spaceId).addEntity(fileInfo, link)` payload encoding — either send `fileInfo.id` (UUID) where the server expects UUID, or accept name-based references and resolve server-side.

### Suggestions for the scenario
- None from this run; auto-unblock probe handles the env transition cleanly.
