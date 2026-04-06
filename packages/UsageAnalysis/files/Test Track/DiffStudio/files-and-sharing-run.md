# Files & Sharing in Diff Studio — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Navigate to Browse > Files > App Data > Diff Studio > Library; click pk.ivp | PASS | PASSED | Navigated to `/files/System.AppData/DiffStudio/library?browse=files`; clicked `pk.ivp` link; PK model opened in DiffStudio (3 cols: t, depot, centr; 1201 rows; step=0.01, count=1) |
| 2 | Set Step to 0.1, Count to 4; view updates | PASS | PASSED | Used native setter + `change` event on step input (0.01→0.1) and count input (1→4); URL updated automatically; charts updated to 4 cyclic dose cycles; rows changed from 1201 to 484 |
| 3 | Copy URL, open new tab, paste URL → same model/inputs load | PASS | PASSED | URL: `.../pk.ivp?params:begin=0&end=12&step=0.1&dose=10000&count=4&...`; new tab loaded same PK model with identical inputs (step=0.10, count=4, 484 rows, 4 cyclic charts) |
| 4 | REMARK: 1-3 curves → only linechart (no Multiaxis & Facet) | PASS | N/A | Confirmed: PK model (2 curves: depot, centr) shows as individual line charts — no Multiaxis/Facet tabs; those appear only for models with 4+ output variables |

## Summary

All 3 functional steps passed. The pk.ivp file opened directly from the Files browser as a DiffStudio view. Modifying Step (0.01→0.1) and Count (1→4) updated the URL query parameters and regenerated the charts immediately. Opening the copied URL in a new tab reproduced the exact same model state. The remark about Multiaxis/Facet absence is confirmed — they only appear for models with sufficient output curves.

## Retrospective

### What worked well
- Direct URL navigation to the library folder (`/files/System.AppData/DiffStudio/library?browse=files`) reliably shows library files
- Clicking `pk.ivp` in the file list opens DiffStudio with the PK model automatically
- URL parameters (`params:step=0.1&count=4`) encode all model inputs — sharing works perfectly by copying the URL
- The new tab loaded exactly the same model state without any login friction (session cookie preserved)

### What did not work
- Slider interaction (as mentioned in the scenario) was not available via browser automation — used native input setter + `change` event dispatch instead (equivalent result)
- `mcp__chrome-devtools__hover` is unreliable for triggering tooltip; dispatching `MouseEvent('mouseover')` works better

### Suggestions for the platform
- The URL parameter format (`params:key=value`) is non-standard; consider documenting this as a deep-link API for DiffStudio
- A "Copy link" button on the DiffStudio ribbon would make the sharing workflow more discoverable

### Suggestions for the scenario
- Step 2 says "Set Step to 0.1 using the slider" — in the actual UI, Step is a text input with h units, not a slider. The scenario should say "text input" or "field"
- Step 3 could specify that parameters are encoded in the URL query string (so users know what to expect before copying)
- Add verification: confirm that the row count in the new tab matches the original tab (484 rows for step=0.1, count=4)
