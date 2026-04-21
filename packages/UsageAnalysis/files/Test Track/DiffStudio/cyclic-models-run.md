# Cyclic Models in Diff Studio — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open Diff Studio + PK-PD Model (Library > PK-PD) | 8s | PASS | PASSED | Invoked `DiffStudio:pkPdNew`; view 'PK-PD' opened; 18 input hosts present (interval, dose, count, central, peripheral, rate-constant, clearance, central-volume, inter-rate, peri-volume, effect, Rate, begin, step, depot, init-effect). |
| 2 | Both Multiaxis and Facet plots updated | 1s | PASS | FAILED | 2b: "Multiaxis" and "Facet" tab text both visible in body. Spec run: body-text scan returned false — likely chart tabs not yet rendered or rendered inside canvas/shadow-unreachable text. Not patched per strict sanity-pass rule. |
| 3 | Modify Count input via clickers; real-time update | 1s | PARTIAL | PASSED | `[name="input-host-count"]` exists. Live `.ui-input-plus`/`.ui-input-minus` click + chart re-render not exercised in 2b; spec only asserts host presence. |
| 4 | Tooltips on Begin, End, Step | 1s | AMBIGUOUS | PASSED | `begin` and `step` hosts present; no `end` host exists on PK-PD. Tooltip hover+wait not exercised. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m |
| grok-browser execution (scenario steps) | 10s |
| Execute via grok-browser (total) | 2m |
| Spec file generation | 30s |
| Spec script execution | 14s |
| **Total scenario run (with model)** | ~3m |

## Summary

Reproduced the Cyclic Models (PK-PD) scenario via `DiffStudio:pkPdNew`. Step 1 passed cleanly; Step 2 passed in the MCP run (tab labels visible) but the standalone Playwright body-text scan failed — the scan likely ran before tabs rendered, or the tab labels live inside canvas/shadow trees not surfaced by `innerText`. Steps 3 and 4 are PARTIAL/AMBIGUOUS by design: live clicker re-render and tooltip hover were not exercised in 2b. **Total scenario run (with model)**: ~3m.

## Retrospective

### What worked well
- Direct `DG.Func.find({package: 'DiffStudio', name: 'pkPdNew'}).apply({})` reliably opens the PK-PD view without traversing the Diff Studio hub / ribbon-combo / Library menu cascade.
- 18 input hosts become queryable via `[name="input-host-*"]` immediately after the view activates.

### What did not work
- Body-text scan for "Multiaxis"/"Facet" in the Playwright run — inconsistent with the MCP observation. Root cause likely timing (chart tabs render asynchronously) or tab labels living in a region not captured by `document.body.innerText`.
- Scenario's "End" input does not exist on PK-PD (only `begin` and `step`), making Step 4 partially unverifiable as written.

### Suggestions for the platform
- Input captions on PK-PD (and similar cyclic models) are lowercase + no-space (e.g. `begin`, `step`, `init-effect`) — consider normalizing to the scenario's human-readable form (`Begin`, `Step`) for consistency across models and easier cross-scenario automation.
- Expose stable `aria-label` or data attributes on Multiaxis/Facet tab elements so body-text scans are not required to detect tab presence.

### Suggestions for the scenario
- Step 4 references **Begin/End/Step** inputs, but the PK-PD model exposes `begin` and `step` (lowercase) with no `end` input — update the scenario to match actual input captions (or request that the platform normalize them; see platform suggestion above).
- Step 3 should explicitly name the clicker controls (`.ui-input-plus` / `.ui-input-minus`) and define an observable outcome (row count / chart range change) so automation can assert the "real-time update" claim.
- Step 1 should mention the direct entry point (`DiffStudio:pkPdNew`) as an automation-friendly alternative to the Apps > Diff Studio > Open model > Library > PK-PD chain.
