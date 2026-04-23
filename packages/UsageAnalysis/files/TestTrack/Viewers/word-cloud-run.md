# Word cloud — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI dataset (setup + semType wait) | 18s | PASS | PASSED | 3624 rows, 88 cols; Bio/Chem wait applied (Molecule columns detected: Structure, Core, R1–R3, R100, R101) |
| 2 | Open Add Viewer gallery, click Word Cloud tile | 4s | PASS | PASSED | Gallery opens with "Add Viewer" title; `.d4-item-card.viewer-gallery` with innerText "Word Cloud"; viewer added (`tv.viewers` includes type "Word cloud"). No viewer errors; only unrelated 404s for missing icon resources |
| 3 | Close viewer, re-add via Toolbox icon | 2s | PASS | PASSED | Close via `[name="Close"]` on `.panel-base`; re-add via `[name="icon-Word-cloud"]` (capital W) in Toolbox Viewers section |
| 4 | Interact with hamburger menu (displayed-text clicks skipped) | 2s | PASS | PASSED | Hamburger `[name="icon-font-icon-menu"]` opens menu with items: General, Pick Up / Apply, Properties..., To Script. Canvas word clicks not automatable from DOM |
| 5 | Open Property Pane via Gear icon | 2s | PASS | PASSED | Gear `[name="icon-font-icon-settings"]` opens Context Panel; Data category visible in `.property-grid-category` |
| 6 | Modify column=Stereo Category, minTextSize=20, maxTextSize=80, bold, rotationStep=90 | 3s | PASS | PASSED | All 5 property setters applied via `wc.props.*`; viewer re-rendered. Valid prop names: `columnColumnName`, `shape`, `minTextSize`, `maxTextSize`, `minRotationDegree`, `maxRotationDegree`, `rotationStep`, `gridSize`, `drawOutOfBound`, `fontFamily`, `bold` |
| 7 | Cleanup: close viewer | 1s | PASS | PASSED | Clean close via `[name="Close"]`; `tv.viewers` no longer contains Word cloud |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 45s |
| grok-browser execution (scenario steps) | 32s |
| Execute via grok-browser (total) | 2m 17s |
| Spec file generation | 50s |
| Spec script execution | 27s |
| **Total scenario run (with model)** | 4m 5s |

## Summary

All 7 scenario steps PASS in both the MCP run against https://dev.datagrok.ai and the
generated Playwright replay (27.2s, all softSteps green). Word Cloud adds via both entry
points (Add-Viewer gallery and Toolbox `[name="icon-Word-cloud"]`), hamburger menu opens
with expected items, gear opens Property Pane, and all five property modifications
propagate. No Word-cloud-specific console errors; only 4 unrelated 404s for missing icon
resources unrelated to the viewer.

## Retrospective

### What worked well
- Add-Viewer gallery click is reliably driven by `.d4-item-card.viewer-gallery` filter on
  `innerText === 'Word Cloud'` rather than relying on the Playwright `text=` operator
  (which can match multiple rows when "Recently used" also contains the viewer)
- `wc.getProperties().map(p => p.name)` from inside `page.evaluate` quickly enumerates the
  valid property surface and avoids guessing names
- The shared `spec-login.ts` + `softStep`/`stepErrors` wrapper kept the spec compact and
  failure-isolating — every softStep recorded independently

### What did not work
- The name `columnColumnName` is still counterintuitive for a single-column viewer — the
  obvious guess is `column` or `wordColumnName`. Reflection-generated naming leaks into
  the public property surface
- Canvas-based word rendering prevents the Step 4 sub-requirement "pressing the displayed
  text should select appropriate rows" from being exercised via DOM automation — we can
  only verify the hamburger menu and submit a screenshot
- Scenario step 1 title reads "World cloud" (typo) in the original `.md` — cosmetic but
  carries through copy-paste into other artifacts

### Suggestions for the platform
- Rename `columnColumnName` to `column` for viewers whose sole input is one column — the
  double-`Column` suffix is an auto-generated artifact that leaks into the public API
- Add a `name=` attribute to the "Add viewer" ribbon icon so selectors don't have to fall
  back to `aria-label` (e.g. `[name="icon-add-viewer"]`)
- Emit a curated list of valid `shape` enum values (currently only `circle` default is
  visible) in the property-grid tooltip/help

### Suggestions for the scenario
- Fix the "World cloud" → "Word cloud" typo in the scenario title and body
- Step 4 "Pressing the displayed text should select appropriate columns" is ambiguous —
  clicking words in a Word Cloud filters rows by that word, it does not change columns.
  Clarify to "…should filter the dataframe to rows matching that word"
- Step 5 mentions the **Property Pane** — explicitly say it appears in the **Context
  Panel (right side)** in Tabs mode, since the Toolbox on the left contains unrelated
  sections
- Step 6 "Modify various properties" is open-ended — enumerate a concrete suggested set
  (e.g. `Column`, `Min text size`, `Max text size`, `Bold`, `Rotation step`) to make the
  test deterministic
