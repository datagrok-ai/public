# Projects / Lifecycle Script: provisioned df-output script source — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 1m 54s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 54s |
| **Total scenario run** | 1m 54s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 54s.

Lifecycle: provision a df-output script → call it → save with `.script` provenance → reopen → re-execute the script.

## Retrospective

### What worked well
- Provisioned-script lifecycle (test fixture script created via `dapi.scripts.save` then called) produces the canonical `Var = ScriptNamespace:ScriptName(...)` provenance and reopens correctly.

### What did not work
- Nothing notable.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- None from this run.
