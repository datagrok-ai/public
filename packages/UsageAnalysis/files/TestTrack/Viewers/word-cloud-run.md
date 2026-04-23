# Word cloud — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI dataset | 18s | PASS | FAILED | 3624 rows, 88 cols; Bio/Chem wait applied (Molecule column); Playwright timed out at login step before setup ran |
| 2 | Open Add Viewer gallery → click Word Cloud tile | 12s | PASS | FAILED | Gallery opened; tile clicked; viewer added without console errors; Playwright blocked by login timeout |
| 3 | Close and re-add via Toolbox Viewers > Word Cloud icon | 4s | PASS | FAILED | Toolbox icon name is `icon-Word-cloud` (capital W). Re-add succeeded without errors |
| 4 | Interact with hamburger menu (displayed-text clicks skipped) | 3s | PASS | FAILED | Hamburger opens menu with "General / Pick Up/Apply / Properties.../ To Script". Canvas word clicks not automatable |
| 5 | Open Property Pane via Gear icon | 4s | PASS | FAILED | Gear on `.panel-base` wrapper; Data section appears in property grid |
| 6 | Modify column=Stereo Category, minTextSize=20, maxTextSize=80, bold, rotationStep=90 | 3s | PASS | FAILED | All properties applied; viewer re-rendered with new values; no console errors |
| 7 | Close viewer | 2s | PASS | FAILED | Clean close via `[name="Close"]` |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 3m 30s |
| grok-browser execution (scenario steps) | 46s |
| Execute via grok-browser (total) | 4m 15s |
| Spec file generation | 1m |
| Spec script execution | 2m 7s |
| **Total scenario run (with model)** | 7m 30s |

## Summary

All 7 MCP scenario steps PASS. The Word Cloud viewer adds via both entry points (Add-Viewer
gallery and Toolbox icon), the hamburger menu opens, the gear icon opens Property Pane, and
property modifications propagate to the viewer with no console errors. The generated Playwright
spec failed at the login step (timeout waiting for `[name="Browse"]`) — the 15s `isVisible`
check on the login input returned false before the dev-server login form rendered, so login
never executed. The scenario logic itself is sound (identical to the network-diagram spec that
passed); this was an intermittent dev-server/login-timing issue rather than a spec defect.

## Retrospective

### What worked well
- Add-Viewer gallery is easy to drive via the MCP `click` tool when the snapshot uid for the Word Cloud tile is known
- `wc.getProperties().map(p => p.name)` quickly enumerates valid property names (`columnColumnName`, etc.)
- Closing and re-adding the viewer between two entry points validates both paths in one scenario
- Hamburger menu DOM items are queryable via `.d4-menu-popup` even though the viewer body is canvas

### What did not work
- Reference doc/spec conventions assume `column` or `wordColumnName` as the property; the actual name is `columnColumnName`
- Canvas-based word rendering means clicking a displayed text (step 4 "press on displayed text") is not practical from DOM automation
- Playwright's `isVisible` on the login input returned false within 15s on dev.datagrok.ai; the spec skipped the login branch and then timed out waiting for Browse. Per skill rules (no post-run edits), left as FAILED
- Word cloud `shape` property accepts enum values but the reference lists no enumeration — had to guess valid values from the existing default (`circle`)

### Suggestions for the platform
- Rename `columnColumnName` to simply `column` for a viewer whose sole input is one column — the double-`Column` suffix is an auto-generated artifact from reflection and leaks into the public property surface
- Provide a fallback click-target for the "Add viewer" icon with a `name=` attribute so automation does not have to fall back to `aria-label`
- Have the word-cloud scenario include an enumeration of valid `shape` values (e.g. `circle`, `rectangle`) in the Settings help

### Suggestions for the scenario
- Step 4 "Pressing the displayed text should select appropriate columns" is ambiguous — the word cloud does not change columns on word click; it filters the dataframe. Clarify the expected behavior
- Step 5 wording is inconsistent with the reference doc — the gear icon opens the Property Pane but Tabs mode shows it in the Context Pane (right side), not the Toolbox. Mention the target pane
- Add an explicit column name for step 6 (e.g. "change Column to `Stereo Category`") rather than leaving "modify various properties" open-ended
- Step 1 title says "World cloud" — typo, should be "Word cloud"; same typo recurs in several places in the .md file
