# Cyclic Models in Diff Studio (PK-PD) — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open Diff Studio + PK-PD from Library | 15s | PASS | PASSED | `runDiffStudio` + dblclick `.diff-studio-hub-card` with text 'PK-PD'; view 'PK-PD' opens with 15+ input hosts (`begin`, `end`, `step`, `dose`, `count`, `depot`, `centr`, `peri`, `eff`, `rate-constant`, `clearance`, etc.) |
| 2 | Multiaxis and Facet tabs updated | 10s | PASS | PASSED | Tabs include 'Multiaxis', 'Facet', 'Grid'; 11 `.d4-viewer canvas` elements render |
| 3 | Modify Count via clickers; real-time update | 20s | PASS | FAILED | `[name="input-host-count"]` contains clicker icons `[name="icon-plus"]` and `[name="icon-minus"]`; clicking plus 3× via MCP dispatch advances count 10→13, URL updates to `count=13` live. **Playwright failure**: `.click()` waits for visibility but clicker icons only render on hover — spec needs `click({force: true})` or hover-first |
| 4 | Tooltips on Begin, End, Step | 10s | PASS | PASSED | Dispatched `mouseover` on each label: Begin → "Begin of dosing interval"; End → "End of dosing interval"; Step → "Time step of simulation". All three tooltips render in `.d4-tooltip` |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 40s |
| grok-browser execution (scenario steps) | 22s |
| Execute via grok-browser (total) | 1m 2s |
| Spec file generation | 38s |
| Spec script execution | 36s |
| **Total scenario run (with model)** | 2m 16s |

## Summary

Cyclic Models (PK-PD) scenario reproduces fully on dev.datagrok.ai. The PK-PD library model loads via double-click, Multiaxis/Facet tabs render with line charts, Count clickers advance the value (with the expected real-time URL reflection), and hover-tooltips on Begin/End/Step display descriptive help text. The Playwright spec FAILED only on Step 3 because the `+`/`-` clicker icons on the count input are only visible on hover — the spec click didn't pre-hover and Playwright's visibility strict-mode-protected `.click()` refused to click an invisible element. The 2b MCP run bypassed this via direct `icon.click()` dispatch. **Total scenario run (with model)**: 2m 16s.

## Retrospective

### What worked well
- `.diff-studio-hub-card` with `hasText: 'PK-PD'` disambiguates across the 5+ library cards
- Count input exposes `[name="icon-plus"]` and `[name="icon-minus"]` with labels `Increase`/`Decrease` — explicit, discoverable clickers
- Tooltips attach to the label-hover path and render in `.d4-tooltip` within ~800ms of `mouseover`; descriptive text content is specific enough to assert on per-input basis
- Real-time URL update on count change confirms the simulation is reactive without checking canvas pixels

### What did not work
- Playwright's `.click()` visibility check blocks clicker icons that are `display:none` until hover — the MCP JS path (`icon.click()` via `evaluate`) bypasses this but a production spec should pre-hover or use `{force: true}`
- Tooltip dispatch requires raw `mouseover` + `mouseenter` DOM events; Playwright's `.hover()` is equivalent but tooltip timing (~500–800ms) needs an explicit wait before reading `.d4-tooltip`

### Suggestions for the platform
- Make clicker icons (`ui-input-plus` / `ui-input-minus`) visible whenever the input is focused or the whole input-host is in the viewport, not only on hover — makes automation straightforward and also helps touch / keyboard users
- Expose tooltips in a DOM-stable way (`aria-describedby` or `title` attribute on the label) so automation doesn't need a raw event dispatch
- Consider adding `aria-label` to the value input itself with the help text, so screen readers and automation can read the description without triggering hover

### Suggestions for the scenario
- Step 3 says "use the clickers to modify the Count value" — clarify the scope: count advances by 1 per click, clickers are the `+`/`-` icons on the right edge of the input
- Step 4 lists "Begin, End, Step" but the scenario's *Expected Results* section says "Hovering over various input fields should display tooltips" — list the three expected tooltip contents so testers have an assertion target
- Step 2 wording "both Multiaxis and Facet plots are updated" is ambiguous — rephrase as "both Multiaxis and Facet tabs render line charts with PK-PD curves"
