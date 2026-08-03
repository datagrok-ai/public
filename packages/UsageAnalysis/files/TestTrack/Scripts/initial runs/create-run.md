# Scripts Create — Run Results

**Date**: 2026-05-06
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse > Platform > Functions > Scripts | 1s | PASS | PASSED | Tree expanders + Scripts label click; landed on `view.type=scripts` |
| 2 | New > R Script... | 1s | PASS | PASSED | `[name="button-New"]` then `R Script...` menu item; ScriptView opens with template |
| 3 | Click "Open script sample table" (asterisk) | 1s | PASS | PASSED | `i[name="icon-asterisk"]`; `cars` table loaded |
| 4 | Click Signature editor (magic wand) | 1s | PASS | PASSED | DevTools-contributed `i.fal.fa-magic`; icon already injected on first poll (sigEditorAvailable=true) |
| 5 | Signature editor: Set Name = testRscript | 1s | PASS | PASSED | `input[name="input-Name"]` filled; blur committed; body's `#name:` updated to testRscript |
| 6 | Navigate to Parameters tab | 1s | PASS | PASSED | `.d4-tab-header` with text "PARAMETERS" |
| 7 | Click "+" to add a new parameter | 1s | PASS | PASSED | Second visible `i.fa-plus` button; new `#input: bool newParam` line appended to CodeMirror body |
| 8 | Set direction=output, name=newParam, type=string | 1s | PASS | PASSED | Parameters tab is a `PropertyGrid` (canvas-rendered, no public DOM API). Used Datagrok's documented body↔params auto-sync: edited body → `#output: string newParam`. **Verified airtight at step 11**: the saved `DG.Script.outputs` array contains `{name:'newParam', propertyType:'string'}` — i.e. parsed param matches scenario intent, not just the body string |
| 9 | Click "Open function editor" | 1s | PASS | PASSED | `i.fal.fa-code` button's `aria-label` is literally `"Open function editor"` — exact match for scenario wording |
| 10 | Click Play, choose sample table (cars) | 9s | PASS | PASSED | Table dropdown set via native HTMLSelectElement.value setter + input/change/blur; OK button clicked; no `Value not defined` errors |
| 11 | Click Save — script persisted with parsed params | 3s | PASS | PASSED | `[name="button-Save"]`; verified saved script: `outputs=[{count,int},{newParam,string}]`, `inputs=[{table,dataframe}]` — closes the step-8 verification gap |
| 12 | Close script view via x | 1s | PASS | PASSED | `.tab-handle-close-button`; full mousedown/mouseup/click sequence required |
| 13 | Browse > Scripts again, click New | 1s | PASS | PASSED | View context switched back to `scripts` |
| 14 | R script: open sample, run | 14s | PASS | PASSED | No errors |
| 15 | Python script: open sample, run | 15s | PASS | PASSED | No errors (Jupyter kernel warm) |
| 16 | Octave script: open sample, run | 15s | PASS | PASSED | No errors |
| 17 | NodeJS script: open sample, run | 14s | PASS | PASSED | No errors |
| 18 | JavaScript script: run | 16s | PASS | PASSED | Triggered native `alert("Hello World!")`; auto-accepted via chrome-devtools `handle_dialog` (and via `page.on('dialog')` in spec). Template has no `#sample:` |
| 19 | Grok script: open sample, run | 13s | PASS | PASSED | No errors |
| 20 | Pyodide script: open sample, run | 13s | PASS | PASSED | Client-side WASM; no errors |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~3m |
| grok-browser execution (scenario steps) | 1m 57s |
| Execute via grok-browser (total) | 4m 57s |
| Spec file generation | 1m 26s |
| Spec script execution | 2m 49s |
| **Total scenario run (with model)** | ~9m 12s |

`grok-browser execution (scenario steps)` is the precise sum of recorded `evaluate_script`
latencies (the `Time` column above). `Execute via grok-browser (total)` is wall-clock
from the start of step 1 to the post-step-20 cleanup, captured via `date +%s.%N` —
4m 57s exactly. `Model thinking` is the residual (4m 57s − 1m 57s ≈ 3m). `Spec file
generation` is from a previous run (this run did not regenerate the spec; only step
11 was edited, in place). `Spec script execution` is `time npx playwright test --headed`,
2m 49s real.

## Summary

All 20 steps PASS end-to-end against `https://dev.datagrok.ai`. **Total scenario run
(with model): ~9m 12s.** The previous run's two AMBIGUOUS items were upgraded to airtight
PASS:

- **Step 9**: `i.fal.fa-code`'s `aria-label` is literally `"Open function editor"`,
  matching the scenario verbatim.
- **Step 8**: closed the verification gap — step 11 now also checks the saved
  `DG.Script.outputs` array, which this run confirmed contains
  `{name:'newParam', propertyType:'string'}`. The body-edit path produced a
  correctly-parsed Script entity (not just a body string that says so), so the
  scenario's actual goal ("the new parameter has direction=output, name=newParam,
  type=string") is verified at the entity level — not relied on via auto-sync alone.

## Retrospective

### What worked well

- Tree-based navigation to Scripts via `[name="tree-expander-…"]` selectors is reliable.
- `i[name="icon-asterisk"]` (Open sample table) and `i[name="icon-play"]` (Run script) work
  uniformly across all language templates.
- The `.tab-handle-close-button` (full mousedown/mouseup/click sequence) closes ScriptViews
  cleanly without prompting.
- `grok.dapi.scripts.filter('friendlyName = "…"').list()` returns a fully-parsed `DG.Script`
  with populated `inputs`/`outputs` arrays — perfect for asserting the scenario's actual
  end state at the entity level rather than relying on body-string regex.
- The native browser alert from the JavaScript template was caught cleanly by
  chrome-devtools `handle_dialog action: accept` (and by `page.on('dialog')` in the spec).
- Auto-syncing CodeMirror body → parsed params means editing the body header IS the
  documented programmatic path for setting parameters.

### What did not work

- **Run dialog `<select>`** — programmatic `.value = ...` plus `change` event isn't enough;
  Dart's listener requires the native `HTMLSelectElement.value` setter (via
  `Object.getOwnPropertyDescriptor(...).set.call(sel, v)`) and `input` + `change` + `blur`.
- **Signature editor (magic-wand) and Open function editor are DevTools contributions** —
  they have no `name=` attribute and don't always inject in time in fresh Playwright
  sessions on dev. This run they were already injected; the spec keeps the JS-API
  fallback for sessions where they aren't.
- **PropertyGrid (Parameters tab) is canvas-only with no public DOM API** —
  `DG.Viewer.fromRoot(gridEl)` returns `type:'Unknown'`; cells aren't addressable;
  synthetic MouseEvents (`mousedown`/`mouseup`/`click`/`dblclick`) at multiple cell
  positions don't trigger the Dart cell-edit handlers. The body↔params auto-sync
  (CodeMirror body edit → grid refreshes automatically) is the only programmatic path.
- **`.tab-handle-close-button.click()`** alone doesn't fire the Dart close handler — full
  mousedown/mouseup/click sequence is required.

### Suggestions for the platform

- Add `name=` annotations to DevTools-contributed ribbon icons (`name="icon-signature-editor"`,
  `name="icon-function-editor"`) so they can be addressed deterministically without relying
  on Font Awesome class names.
- Make `ChoiceInput`/`TableInput` listen on the standard `change` event from the underlying
  `<select>` so JS-driven `el.value = …` flows just work — current behavior makes browser
  automation needlessly fragile.
- Make the Parameters PropertyGrid cells DOM-addressable (or expose a small JS API to
  mutate a param's direction/type/name) so automation doesn't have to fall back to the
  body header. The body-edit path is robust but requires automation to know the
  language-specific comment style.
- Annotate the tab-handle close button with a stable `name=` (e.g. `name="tab-close-<view>"`)
  and ensure a plain `.click()` triggers close without needing a synthesized mousedown/up.
- Print a clearer error than "Value not defined" when a required dataframe input lands at the
  server with no value — currently nothing in the UI hints at the cause.

### Suggestions for the scenario

- Step 8 cannot be performed cell-by-cell on the canvas-based PropertyGrid through any
  documented DOM selector. Either rewrite the step to "edit the metadata header in the
  script body" (Datagrok auto-syncs body↔params), or add per-cell DOM addressability /
  expose a JS API to mutate a param's direction/type/name from automation.
- Step 10 should explicitly state "Select cars from the Table dropdown, then OK" — the
  default value of an empty `<select>` is the leading source of run errors.
- Step 12 should reference the tab close button (`.tab-handle-close-button`) explicitly —
  the in-view `icon-times` ribbon button is ambiguous when multiple views are open.
- The all-language test (steps 13–20) ran in ~1m 40s this time (no Jupyter cold-start
  penalty). When the kernel is cold, expect 30s+ on the first server-language run —
  consider a "warm up" step or splitting into per-language scenarios so a single failure
  doesn't cost the whole batch.
