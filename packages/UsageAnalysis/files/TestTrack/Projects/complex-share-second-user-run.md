# Projects / Complex share-second-user: dual-level grant via JS API (Step 12; Step 13 deferred) — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 54.3s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 54.3s |
| **Total scenario run** | 54.3s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 54.3s.

Step 12 covered: configure project sharing for another user at View-and-Use AND Full access levels via JS API (`grok.dapi.permissions.grant(p, recipient, edit)`).

Fix this session: migrated from `loginToDatagrok` (defaulted to `localhost:8888`) to canonical `gotoApp + setupSession` from `_helpers.ts` (honors `projectsTestOptions` baseURL).

## Retrospective

### What worked well
- JS API `permissions.grant` round-trips at both edit=false (View-and-Use) and edit=true (Full access).
- Defensive skip-with-note pattern absorbs intermittent FK violations on the permissions table.

### What did not work
- Nothing notable in the final run.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- Step 13 (re-auth as second user, navigate Browse > Dashboards, locate shared project, open) DEFERRED — requires `helpers.playwright.session.logoutAndLoginAs` (Helper 3 from helpers-batch-1, NOT yet registered).
