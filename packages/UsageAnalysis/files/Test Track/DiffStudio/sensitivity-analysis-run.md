# Sensitivity Analysis (Diff Studio) — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open Diff Studio + Bioreactor; Edit off; form with inputs opens | 6s | PASS | FAILED | 2b: opened Bioreactor via `DG.Func DiffStudio:demoBioreactor`, view name 'Bioreactor', 23 inputs visible. Spec asserts `[name^="input-host-"]` count > 0 — returned 0 on dev. The Bioreactor model form inputs do not expose `input-host-*` `name=` attributes, only CSS classes. Edit toggle state is not verifiable via any stable selector. |
| 2 | Click Sensitivity icon on ribbon | 5s | PASS | PASSED | Clicked `.diff-studio-ribbon-sa-icon`; new view "Bioreactor - comparison" opened with 30 inputs including Method, Samples, and min/max ranges for each model parameter. |
| 3 | Select Parameters (Process mode, FFox, KKox, FFred) | 3s | PARTIAL | PASSED | Verified `input-host-FFox`, `input-host-KKox`, `input-host-FFred` are present on the SA form. Switcher toggles (Dart-side widgets) not scripted — 2b confirmed the form is populated but did not enable the switches. |
| 4 | Click Run icon on ribbon; 4 viewers open | 60s | FAIL | FAILED | 2b: `.grok-icon.fal.fa-play` located by class; `.click()` did NOT start the analysis — polled `tv.viewers` for 60s, still only `[Grid]`. Spec used `click({force: true})` and polled 90s, still got 1 viewer. Either the click does not reach the Dart handler (needs a real pointer-event sequence) or the SA run requires an enabled switcher first. No balloon, no warning, no error was emitted. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 30s |
| grok-browser execution (scenario steps) | 1m 15s |
| Execute via grok-browser (total) | 3m 30s |
| Spec file generation | 1m |
| Spec script execution | 1m 42s |
| **Total scenario run (with model)** | 6m 12s |

All rows are full-phase wall-clock (incl. model thinking and retries), not just tool latency. The two `scenario steps` rows roughly sum to `Execute via grok-browser (total)`.

## Summary

2 of 4 steps passed cleanly in 2b (steps 1 and 2), 1 was PARTIAL (step 3 — inputs are present but switchers weren't toggled), and 1 failed outright (step 4 — Run click never produced the 4 viewers, after 60s in 2b and 90s in the Playwright spec). The Playwright replay reproduced the failure faithfully: the Sensitivity icon opens the comparison view and the parameter input hosts exist, but clicking `.grok-icon.fal.fa-play` silently no-ops — no balloon, no warning, no error. Step 1 failed under Playwright for a different reason: the Bioreactor model form inputs don't expose `input-host-*` `name=` attributes, so the spec's count assertion returned 0 even though the form is clearly rendered. **Total scenario run (with model)**: 6m 12s.

## Retrospective

### What worked well
- `DG.Func.find({package: 'DiffStudio', name: 'demoBioreactor'})` reliably loads the Bioreactor model.
- `.diff-studio-ribbon-sa-icon` is a stable CSS hook for the Sensitivity icon; the comparison view opens predictably.
- The SA form surfaces `input-host-FFox`, `input-host-KKox`, `input-host-FFred` so parameter presence is testable via DOM.
- Standalone Playwright spec (explicit login, 120s login wait, Tabs mode, 2s settle after `page.goto`) worked end-to-end against dev without any MCP coupling.

### What did not work
- `.grok-icon.fal.fa-play` click does not trigger the SA run. Both `.click()` (2b) and Playwright `click({force: true})` (2e) failed to produce the 4 viewers.
- Bioreactor model form inputs don't use `input-host-*` `name=` attributes (only CSS classes), so DOM-based input enumeration via `[name^="input-host-"]` returns 0 on that view — a selector mismatch rather than a missing form.
- Parameter switchers for FFox / KKox / FFred on the SA form were not scriptable from 2b; toggling them likely requires exercising the Dart-side widget (real pointer events) rather than a synthetic DOM event.
- The Run button provides no user feedback when it no-ops — indistinguishable from a hung analysis.

### Suggestions for the platform
- Diff Studio ribbon icons (Sensitivity, Fit, Run play) need stable `name=` / `aria-label` selectors. Today they are distinguishable only by CSS class and FontAwesome glyph, which is fragile for automation and accessibility.
- Run icon click on the SA ribbon may silently no-op if no parameter switcher is enabled — emit a balloon/warning when Run is pressed without any selected parameters, so users (and automation) get actionable feedback.
- Expose `name="input-host-<param>"` attributes uniformly on both the Bioreactor model form and the SA form; right now the main model view uses CSS-only labels while the SA form already has proper `name=` hosts.
- Provide a scriptable toggle for SA parameter switchers (e.g. `ui-input-switch` with a data-param attribute, or a JS API shortcut like `sa.enableParams(['FFox','KKox','FFred'])`).

### Suggestions for the scenario
- Step 1 "Edit toggle off" is not observable from the DOM today — reword to describe the expected initial form state or drop the toggle check.
- Step 3 should spell out that enabling a parameter requires flipping its switcher — and that min/max range inputs are inert until the switcher is on.
- Step 4 should clarify the prerequisite ("at least one parameter switcher enabled") and name the expected four viewers, so a tester knows what to look for when the Run click appears to no-op.
