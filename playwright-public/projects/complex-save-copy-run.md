# Projects / Complex Save Copy: round-trip Save Copy with sync OFF/ON — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 3m 06s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 3m 06s |
| **Total scenario run** | 3m 06s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 3m 06s.

Steps 1-7 covered: open demog from File share; Save baseline with Sync ON; Save Copy with Sync OFF via inline UI mode-radio + Data Sync toggle; Close + reopen NoSync copy; Re-save under `<root>_Sync` with Sync ON; Verify post-resave state via reopen + row-count parity vs underlying CSV; 3-project best-effort cleanup in finally.

Step 6 "force Data Sync refresh" surface reduced to reopen-driven re-materialization + row-count parity (UI selector for per-table refresh trigger undocumented in references — see prior SCOPE_REDUCTION proposal #1).

## Retrospective

### What worked well
- Save Copy mode-radio chooser drives `'original' | 'copy' | 'personal-view-customizations'` correctly via label-click + radio-input synchronization.
- 3-project finally-cleanup leaves no orphans on dev.

### What did not work
- Nothing notable.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- Step 6 per-table Data Sync refresh trigger UI selector still undocumented; current spec verifies via reopen + row-count parity.
