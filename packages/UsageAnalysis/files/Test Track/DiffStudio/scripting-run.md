# DiffStudio Scripting — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open Diff Studio, Edit toggle, click </> to open JS script view | 40s | AMBIGUOUS | PASSED | Bioreactor opened via `DiffStudio:demoBioreactor`, but the ribbon "Edit" toggle and "</>" JS view icon render as `div.image-icon.diff-studio-svg-icon` with only a `background-image` URL — no `name`, `aria-label`, `data-*`, or `title`. No reliable DOM selector available. Spec records the blocker via `test.info().annotations` and asserts ribbon icons count > 0 only. |
| 2 | Run the script; adjust "Final at" slider | 20s | AMBIGUOUS | PASSED | Blocked by Step 1 — no ScriptView to run or find `input-Final` in. Spec probes for `ScriptView`, `.fa-play`, and `input[name="input-Final"]` and logs all three as absent. |
| 3 | Add `//tags: model`; Save the script | 15s | AMBIGUOUS | PASSED | Blocked by Step 1 — CodeMirror instance is not mounted without the JS script view. Spec logs `.CodeMirror` absence. |
| 4 | Access Model in Model Hub (Apps > Run Model Hub) | 25s | AMBIGUOUS | PASSED | No `ModelHub` / `modelHub` entry in `DG.Func.find({tags:[FUNC_TYPES.APP]})`. Scenario says "Apps > Run Model Hub" but there is no canonical app function by that name on dev. `Compute2:modelCatalog` is a plausible entry point but disagrees with the scenario wording. |
| 5 | Interact with Model in Model Hub (adjust Final) | 20s | AMBIGUOUS | PASSED | Blocked by Step 4 — cannot open a model from a hub that couldn't be opened by name. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m |
| grok-browser execution (scenario steps) | 4s |
| Execute via grok-browser (total) | 2m 30s |
| Spec file generation | 30s |
| Spec script execution | 15s |
| **Total scenario run (with model)** | 3m 15s |

## Summary

Scenario FAILs on the platform side: every scripted step after "open Bioreactor" is blocked
by UI that is not automation-addressable. The DiffStudio ribbon uses unnamed CSS-only icons
and there is no discrete Model Hub app function to target. The Playwright spec technically
passes (PASSED) because it only encodes presence probes and annotates each blocker as
AMBIGUOUS — it is kept as a documentation target and regression sentinel for future
platform fixes. Total scenario run (with model): 3m 15s.

## Retrospective

### What worked well
- Opening the Bioreactor model directly via `DG.Func.find({package:'DiffStudio', name:'demoBioreactor'})` is reliable and fast.
- `grok.shell.v.name === 'Bioreactor'` is a stable view-identity check.

### What did not work
- DiffStudio ribbon icons cannot be selected — `div.image-icon.diff-studio-svg-icon` carries only a CSS `background-image` URL; no `name`, `aria-label`, `data-*`, or `title`. Root cause: ribbon icons are painted as CSS sprites rather than labeled widgets.
- Model Hub has no canonical app function — `DG.Func.find({tags:['app']})` on dev does not expose a `ModelHub` entry, so scenarios cannot reference it by name.
- Downstream steps (Run, Save, Final slider in hub) cascade-fail once the Edit toggle and `</>` icon cannot be reached.

### Suggestions for the platform
- Add `name=` attributes or proper ARIA roles to DiffStudio ribbon icons (Edit toggle, `</>` JS view, Save to Model Hub, Sensitivity, Fit). Right now icons are distinguished only by `background-image` URL, which is fragile and non-automation-friendly.
- Expose Model Hub as a discrete app function — e.g. `DG.Func.find({tags:['app'], package:'Compute2'})` with a canonical `name: 'ModelHub'` — so scenarios can reference it by name rather than via a menu path.
- Consider `data-diff-studio-action="edit|code|run|save|sensitivity|fit"` on ribbon buttons so automation has a single stable hook.

### Suggestions for the scenario
- Replace obscure ribbon glyph references ("</>", "Edit toggle") with concrete icon `name=` values once the platform adds them, or include screenshots.
- Spell out the exact entry point for Model Hub instead of "Apps > Run Model Hub" — the wording does not match any registered function on dev.
- Explicitly state the pre-condition: "Bioreactor opened via Library" vs "opened via demoBioreactor" — the two paths produce different DOM wrapping.
