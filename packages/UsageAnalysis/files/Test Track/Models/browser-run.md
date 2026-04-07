# Browser — Run Results

**Date**: 2026-03-12
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-----------|-------|
| 1 | Go to Browse > Platform > Predictive models | PASS | PASSED | Navigated to /models. 1 model visible (TestDemog). |
| 2 | Type TestDemog in the search field | PASS | PASSED | Search box functional; TestDemog found. |
| 3 | Check all Context Panel tabs for the model | PASS | PASSED | Tabs visible: Details, Performance, Activity, Sharing, Chats, Sticky meta. Details shows: Created by, Trained on dataset ID, Inputs, Outputs, Method. |
| 4 | Click Filter templates icon and check its content | PASS | PASSED | Filter panel opened; Misc and Description sections visible. |
| 5 | Select multiple models by CTRL+clicking | SKIP | SKIP | Only 1 model existed in the browser; multi-select requires ≥2 models. Could not be tested. |
| 6 | On Context Pane, open Commands tab and click Compare | SKIP | SKIP | Multi-select prerequisite failed; Compare step skipped. |

## Summary

4 of 6 steps passed; 2 skipped because only 1 model was available (TestDemog was the only model — the second numeric model from Train.md was not separately saved). The browser navigation, search, Context Panel tabs, and Filter templates all work correctly. Multi-select and Compare could not be tested.

## Retrospective

### What worked well
- Browse > Platform > Predictive models navigation works correctly
- Search field filters models by name
- Context Panel shows all expected tabs (Details, Performance, Activity, Sharing)
- Filter templates panel opens and shows content

### What did not work
- **Multi-select test required ≥2 models** — the scenario depends on having multiple models available. Only TestDemog was present at this point. The scenario should either (a) train 2 models before Browser.md, or (b) explicitly state that the second model from the "numeric by numeric" training in Train.md must be saved with a name.

### Suggestions for the platform
- Add a keyboard shortcut or tooltip hint for CTRL+click multi-select in the browser grid
- The Compare feature should be documented in the Context Panel Commands tab

### Suggestions for the scenario
- Steps 5-6 require ≥2 models; add a prerequisite note: "Save a second model from the Train scenario before testing multi-select"
- Step numbering skips 5 and 6 (goes 4 → 7 → 8) — fix the numbering
- Add expected descriptions for each Context Panel tab (what should appear in Details, Performance, etc.)
