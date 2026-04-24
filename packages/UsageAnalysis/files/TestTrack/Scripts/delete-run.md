# Scripts Delete — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Platform > Functions > Scripts | 3s | PASS | PASSED | `grok.shell.route('/scripts')` |
| 2 | Find the script via search | 3s | PASS | PASSED | testRscript visible in gallery. In spec: seeded a fresh testRscript via `DG.Script.create(...) + grok.dapi.scripts.save()` because the MCP run had already deleted it |
| 3 | Right-click, select Delete | 5s | PASS | FAILED | contextmenu on `.grok-gallery-grid-item`, `Delete` from menu. Playwright failure: the freshly-saved script was not yet visible in the gallery — `find(...) === undefined` → `getBoundingClientRect` threw on the card. The Scripts gallery does not auto-refresh after a server-side `scripts.save()` |
| 4 | Click YES in confirmation dialog | 3s | PASS | PASSED | Dialog text "Are you sure? Delete script 'testRscript'?". Confirmed via YES button (no `name=` attribute — found by text match) |
| 5 | Verify script is gone | 2s | PASS | FAILED | MCP: `grok.dapi.scripts.filter('name = "testRscript"').list()` returned 0. Playwright failure cascaded — step 3 never clicked Delete, so the seeded script was still present (`remaining: 1`). After the spec, the MCP session cleaned up the orphan seed |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~30s |
| grok-browser execution (scenario steps) | ~25s |
| Execute via grok-browser (total) | 57s |
| Spec file generation | ~25s |
| Spec script execution | 22s |
| **Total scenario run (with model)** | ~1m 45s |

## Summary

Delete passed end-to-end in the MCP session: confirmation dialog was correctly shown ("Are you
sure? Delete script 'testRscript'?"), YES dismissed it, and the script was removed both from
the server (`grok.dapi.scripts.filter(...).list()` returns 0) and the gallery. The Playwright
replay independently seeded a testRscript at step 2 (because the MCP session had already
deleted it), but the gallery did not re-render to include the new card in time, so the
contextmenu step failed. After the spec, the MCP session cleaned up the leftover seeded script
for the next scenario.

## Retrospective

### What worked well
- Context menu on the script card exposes a clear `Delete` entry
- Confirmation dialog message is explicit and names the script
- YES/CANCEL dialog pattern is consistent with other entity deletion flows
- Server-side verification via `grok.dapi.scripts.filter(...).list()` is decisive

### What did not work
- Gallery does NOT refresh after a programmatic `grok.dapi.scripts.save()` — newly saved scripts are invisible until a manual reload or navigation
- Confirmation dialog buttons have no `name=` attribute, so selectors fall back to text match (`YES` / `CANCEL`)
- Search input still doesn't filter the gallery (consistent with browser.md observations)

### Suggestions for the platform
- Scripts gallery should subscribe to `grok.dapi.scripts` change events and re-render when a script is created/deleted server-side
- Add `name="button-YES"` / `name="button-NO"` to confirmation dialogs for deterministic selector targeting
- Consider a Delete keyboard shortcut on the selected card (e.g. `Del`)

### Suggestions for the scenario
- Step 2 "find the script using search" — either describe current behavior (search doesn't filter; relies on visual scanning) or note it as a known gap
- The scenario has no cleanup step; after delete is complete, it should confirm that re-creating the same-named script works (idempotency)
- Specify that the name column is sometimes displayed title-cased (`TestRscript`) but the internal name is the user's original input (`testRscript`)

## Re-run after spec fixes (2026-04-24)

After patching the spec for robust waits (`waitForFunction` on `grok.shell.v?.name`, full
route round-trips to force gallery refresh, Playwright right-click for context menus, JS-API
fallbacks for the Run-dialog table dropdown and the signature-editor's internal state), the
Playwright run now **PASSES** in 17s for Scripts Delete. All scenario steps above that were
previously marked `FAILED` in the Playwright column now pass on the updated spec. Steps still
marked `SKIPPED` are intentional (manual file picker, canvas toolbox, cross-cutting project
flow) and use `test.step.skip` in the spec.
