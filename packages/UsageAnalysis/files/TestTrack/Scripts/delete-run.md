# Scripts Delete — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Platform > Functions > Scripts | 5s | PASS | PASSED | `grok.shell.route('/scripts')` round-trip; `grok.shell.v.name === 'Scripts'` |
| 2 | Find the script via search | 6s | PASS | PASSED | testRscript was seeded via `DG.Script.create + grok.dapi.scripts.save` because the gallery did not contain it. Search input `input[placeholder="Search scripts by name or by #tags"]` accepts the text but does NOT filter the gallery — card was found by polling `.grok-gallery-grid-item-title` |
| 3 | Right-click, select Delete | 5s | PASS | PASSED | Dispatched `MouseEvent('contextmenu', {button:2, clientX, clientY})` on `.grok-gallery-grid-item`; `Delete` from `.d4-menu-item`. Playwright used real `page.locator(...).click({button:'right', force:true})` |
| 4 | Click YES in confirmation dialog | 25s | PASS | PASSED | Dialog text "Are you sure? Delete script 'testRscript'?". Confirmed via `[name="button-YES"]`. **Two iterations needed** — the prior visibility check (`d.offsetParent !== null`) is broken because `.d4-dialog` is `position: fixed`, so `offsetParent` is always null. Fixed visibility test now uses `getBoundingClientRect()` + `getComputedStyle()` |
| 5 | Verify script is gone | 4s | PASS | PASSED | `grok.dapi.scripts.filter('name = "testRscript"').list()` returned 0; gallery (after route round-trip) no longer shows the card |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~30s |
| grok-browser execution (scenario steps) | 15s |
| Execute via grok-browser (total) | 45s |
| Spec file generation | ~30s |
| Spec script execution | 13s |
| **Total scenario run (with model)** | ~1m 30s |

## Summary

All 5 scenario steps PASSED end-to-end in MCP, and the regenerated Playwright spec PASSED on the first attempt in 12.6s. The key change versus the previous run: step 4's dialog-visibility check no longer relies on `offsetParent !== null`, which is always false for fixed-positioned `.d4-dialog` elements and was masking the UI path's success behind a JS API fallback. With the fix, the spec now verifies the actual UI path (right-click → Delete → YES) instead of silently bypassing it.

## Retrospective

### What worked well

- `grok.shell.route('/scripts')` is a clean way to land on the gallery and verify via `grok.shell.v.name === 'Scripts'`
- `[name="button-YES"]` is in fact present on the confirmation dialog (the prior run reported "no name= attribute" — that was wrong; the buttons are annotated)
- Server-side verification via `grok.dapi.scripts.filter(...).list()` is decisive and avoids polling the gallery
- `MouseEvent('contextmenu', {button: 2, clientX, clientY})` reliably opens the Delete menu when the card is in viewport

### What did not work

- Generic visibility helpers using `offsetParent` are unreliable for dialogs/popups with `position: fixed` — `.d4-dialog` matches even when on screen. This is a well-hidden trap because the rest of the DOM behaves the way you'd expect
- The Scripts gallery does not auto-refresh after a server-side `grok.dapi.scripts.save()` — newly saved scripts are invisible until a manual reload or `grok.shell.route` round-trip
- The search input `input[placeholder="Search scripts by name or by #tags"]` accepts text but does not filter the gallery — the card had to be located by scanning all 50 cards

### Suggestions for the platform

- Scripts gallery should subscribe to `grok.dapi.scripts` change events and re-render when a script is created/deleted server-side
- The search input on the Scripts gallery should actually filter visible cards (matching the placeholder text "Search scripts by name or by #tags")
- Consider a Delete keyboard shortcut on the selected card (e.g. `Del`)
- Document somewhere prominent that `.d4-dialog` is `position: fixed` and that selectors / visibility checks must use `getBoundingClientRect()` rather than `offsetParent`

### Suggestions for the scenario

- Step 2 ("find the script using search") — either describe current behavior (search doesn't filter; card has to be scrolled to) or note it as a known gap
- The scenario has no cleanup step; after delete is complete, it could verify that re-creating the same-named script works (idempotency)
- Consider adding an explicit precondition: "If `testRscript` is missing from the gallery, run the Create scenario first (or skip)"
