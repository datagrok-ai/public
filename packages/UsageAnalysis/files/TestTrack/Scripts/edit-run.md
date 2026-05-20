# Scripts Edit — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Platform > Functions > Scripts | 9s | PASS | PASSED | `grok.shell.route('/scripts')` after a `route('/')` round-trip; polled `grok.shell.v?.name === 'Scripts'` and `.grok-gallery-grid-item-title` count > 0 (50 cards) |
| 2 | Find testRscript and double-click | 7s | PASS | PASSED | `.grok-gallery-grid-item-title` text=`testRscript` → closest `.grok-gallery-grid-item` → `dblclick`; polled `grok.shell.v?.type === 'ScriptView'` + `.CodeMirror?.CodeMirror` |
| 3 | Add newParam="test" to script body | 5s | PASS | PASSED | `cm.setValue(before + '\nnewParam="test"\n')` on the active CodeMirror — preserves existing body |
| 4 | Click Save | 8s | PASS | PASSED | `[name="button-Save"]` clicked; persistence confirmed via `grok.dapi.scripts.filter('name = "testRscript"').first()` body containing `newParam="test"` |
| 5 | Close script view | 6s | PASS | PASSED | `.d4-tab-header.selected i[name="icon-times"]` did not match — fell back to `grok.shell.v.close()`; view returned to `Scripts` |
| 6 | Double-click testRscript again; verify newParam="test" | 9s | PASS | PASSED | Re-found gallery card by label, dblclick, polled until `grok.shell.v?.name === 'testRscript'`; CodeMirror body contains `newParam="test"` |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~30s |
| grok-browser execution (scenario steps) | ~14s |
| Execute via grok-browser (total) | 1m 13s |
| Spec file generation | ~10s (no edit needed — existing spec faithfully matches the run log) |
| Spec script execution | 22s |
| **Total scenario run (with model)** | ~1m 45s |

## Summary

End-to-end edit works on dev: the seeded `testRscript` opens, `newParam="test"`
is appended via CodeMirror, `Save` persists the change to the server, and
re-opening the script confirms the new line is present. The Playwright replay
of `edit-spec.ts` passed on the first attempt in 19.6s — every scenario step
green. The only caveat is a missing prerequisite: `testRscript` did not exist
on dev, so the spec's `[pre]` softStep seeded it before step 1 (and we did the
same in the MCP run via `DG.Script.create`).

## Retrospective

### What worked well
- Token exchange `POST /api/users/login/dev/<key>` → `DATAGROK_AUTH_TOKEN` →
  `loginToDatagrok(page)` was clean and matches what `grok test` does internally.
- The existing `edit-spec.ts` was a faithful transcription of the working flow —
  no edits required, ran green on the first invocation.
- `grok.dapi.scripts.filter('name = "testRscript"').first()` is a reliable Save
  verifier (poll the persisted body rather than UI state).
- `.grok-gallery-grid-item-title` text match + dblclick on the closest card is a
  robust alternative to gallery-search filtering (cards have no `name=` attribute).

### What did not work
- Step 5's literal "x icon" selector `.d4-tab-header.selected i[name="icon-times"]`
  did not match — the active tab in the view selector uses different markup. The
  spec works because it skips that path entirely and calls `grok.shell.v.close()`.

### Suggestions for the platform
- Add a stable selector (`[name="view-tab-close"]` or similar) to the active
  view-tab close icon so Test Track scenarios can drive the literal "click x"
  step without needing a JS-API fallback.
- `grok.shell.route()` should expose a promise that resolves when the target view
  is mounted — would replace the `for (let i=0; i<60; i++)` polling loops scattered
  across this and other Scripts specs.
- Browse tree / scripts gallery do not auto-refresh after a JS-API
  `grok.dapi.scripts.save(...)`. The `[pre]` step works around this with a
  `route('/')` → `route('/scripts')` round-trip.

### Suggestions for the scenario
- Add an explicit precondition line: "This scenario assumes `testRscript` already
  exists (created via `Scripts/create.md`). If it was deleted by a previous run
  or this is a fresh server, run create.md first."
- Step 5 says "Click on the **x** icon to close script view" — consider clarifying
  that this means the close icon on the active view tab (not the View → Close
  menu, not Esc), and ideally tagging that icon with a stable `name=` attribute
  so the wording matches the selector.
