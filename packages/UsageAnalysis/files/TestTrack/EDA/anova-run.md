# ANOVA — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog.csv from Demo Files | 1s | PASS | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` + `addTableView`; 5850 rows, 11 cols; semType wait returned fast |
| 2 | Top Menu > ML > Analyze > ANOVA... | 2s | PASS | PASSED | `document.querySelector('[name="div-ML---Analyze---ANOVA..."]').click()`; dialog title "ANOVA" opened; default inputs Category=RACE, Feature=AGE, Alpha=0.05 |
| 3 | Click RUN; Box plot + Analysis + F-test tabs appear | 1s | PASS | PASSED | Button name is `[name="button-Run"]` (mixed case — NOT `button-RUN`/`button-OK`); dialog closes; `grok.shell.tv.viewers` gains Box plot; `d4-tab-host` shows "Analysis" and "F-test" tabs; no console errors, no balloons |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 45s |
| grok-browser execution (scenario steps) | 6s |
| Execute via grok-browser (total) | 1m 30s |
| Spec file generation | 30s |
| Spec script execution | 26s |
| **Total scenario run (with model)** | 2m 26s |

## Summary

All 3 scenario steps passed against dev. Dataset opens via JS API in ~1s; ANOVA dialog
mounts with sensible defaults (RACE, AGE, 0.05); `button-Run` closes the dialog and
produces the expected artifacts — a Box plot viewer attached to the table view and a
`d4-tab-host` with "Analysis" and "F-test" tabs. Expected result from the scenario matches
exactly: "The box plot, Analysis, and F-test tabs are added to the results. There are
no errors in the console." Playwright spec passed in 8.9s (total 26s including browser
launch). **Total scenario run (with model): 2m 26s.**

## Retrospective

### What worked well
- `[name="div-ML---Analyze---ANOVA..."]` click opens the dialog deterministically; a
  single 1.2s wait is enough for layout to settle.
- Defaults on `demog.csv` (Category=RACE categorical, Feature=AGE numeric, Alpha=0.05)
  let Step 3 trigger RUN without any prior input editing.
- Box plot + tab host appearance are observable cheaply — `grok.shell.tv.viewers` enum
  for the box plot and `.d4-tab-host .d4-tab-header` text for the tab labels. Both flip
  well within the 30s soft deadline.
- Scenario text exactly matches observed behavior, including the "no errors in the
  console" clause.

### What did not work
- Nothing — all 3 steps PASS.

### Suggestions for the platform
- The ANOVA dialog's action button is registered as `[name="button-Run"]` (mixed case),
  while most other ML dialogs use `[name="button-RUN"]` or `[name="button-OK"]`.
  Inconsistency in casing of `[name="button-..."]` across ML dialogs (Run vs RUN vs OK)
  is a small but recurring automation footgun — it's the kind of thing that only surfaces
  as a flaky spec step. Standardizing all primary dialog actions on `button-OK` (with
  "Run" as the label) would let spec authors reuse one selector across the whole ML
  dialog family.
- Consider naming the tab host container something like `[name="anova-results"]` or
  `[name="viewer-ANOVA"]` so specs don't have to fall back to scoping via
  `.d4-tab-host` + tab text.

### Suggestions for the scenario
- Step 3 could name the expected viewer explicitly ("Box plot viewer is added to the
  table view") and the tab names verbatim so that spec authors don't need to infer
  them from the screenshot.
- Optional: note the default Category/Feature/Alpha values so readers know no input
  edits are required before clicking Run.

### Suggestions for the skill
- Playwright's default `testMatch` is `*.spec.ts` (dot), but the Test Track convention
  is `*-spec.ts` (dash). The `npx playwright test <path>` invocation silently reports
  "No tests found" without a `testMatch` override. Either adopt `.spec.ts` naming, or
  include a one-off `--config` override pointing at a tiny config file with
  `testMatch: /-spec\.ts$/` in the skill's 2e instructions.
