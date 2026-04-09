# Verify Scripting and Model Interaction in Diff Studio — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open DiffStudio, Edit toggle, click </> icon | PASS | 12s | PASSED | Called `DiffStudio:runDiffStudio`, loaded Bioreactor, clicked Edit span then `</>` span; ScriptView opened with JS code |
| 2 | Run script, adjust Final slider | PASS | 12s | PASSED | Clicked play button via `page.mouse.click()`; js-view-base opened with inputs + table + line chart; Final changed 1000→500 |
| 3 | Add `//tags: model` and Save | PASS | 4s | PASSED | Inserted tag via CodeMirror API; clicked SAVE button; no errors |
| 4 | Access Model Hub, find Bioreactor | PASS | 6s | PASSED | Called `Compute2:modelCatalog`; Bioreactor card visible in catalog |
| 5 | Open model from Hub, adjust Final slider | PASS | 10s | PASSED | Clicked Bioreactor label, called `obj.apply()`; Final changed 500→600; line chart and canvases present |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 50s |
| Spec file generation | 3s |
| Spec script execution | 49s |

## Summary

All 5 steps passed. The DiffStudio Edit toggle and `</>` export icon work correctly, producing a JS script with all model equations and annotations. The script runs and produces a RichFunctionView with inputs, table, and Multiaxis line chart. Adding `//tags: model` and saving makes the script appear in Model Hub. The saved model opens from Model Hub and responds to input changes.

## Retrospective

### What worked well
- Edit toggle clickable via `span.diff-studio-ribbon-text` text match
- `</>` icon found via `span.d4-ribbon-name` with text `</>`
- CodeMirror API (`cm.CodeMirror.getValue()`/`setValue()`) works for programmatic script editing
- `obj.apply()` on the selected Model Hub card opens the model's RichFunctionView
- Inputs use `name="input-{param}"` attributes consistently across DiffStudio and RichFunctionView

### What did not work
- The script's RichFunctionView shows the table but row count verification was unclear — `grok.shell.tables` lists DiffStudio tables, not the script view's internal table
- The REMARK notes "no Process mode" and "just Multiaxis" — confirmed in the RichFunctionView

### Suggestions for the platform
- The `</>` export icon could have an aria-label for accessibility
- The script view could show a confirmation toast/balloon on successful save
- Model Hub could provide a search/filter to find recently saved models quickly

### Suggestions for the scenario
- Step 2 says "Final at" but the input is just "Final" — clarify the name
- Step 3 could mention where `//tags: model` should be inserted (after `//name:` line)
- Step 4 could note the `Compute2:modelCatalog` function for programmatic access
