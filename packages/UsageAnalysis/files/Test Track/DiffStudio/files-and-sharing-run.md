# Files & Sharing in Diff Studio — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Navigate to Browse > Files > App Data > Diff Studio > Library; click pk.ivp | PASS | 25s | PASSED | Navigated to `/files/system.appdata/DiffStudio/library`; clicked `pk.ivp`; PK model opened in DiffStudio (3 cols: t, depot, centr; 1201 rows; step=0.01, count=1) |
| 2 | Set Step to 0.1, Count to 4; view updates | PASS | 8s | PASSED | Typed 0.1 into step field, 4 into count field via UI (Ctrl+A, type, Enter); URL updated; rows changed from 1201 to 484; line chart updated with 4 dose cycles |
| 3 | Copy URL, open new tab, paste URL → same model/inputs load | PASS | 10s | PASSED | URL: `.../pk.ivp?params:begin=0&end=12&step=0.1&dose=10000&count=4&...`; new tab loaded same PK model with identical inputs (step=0.10, count=4, 484 rows) |
| 4 | REMARK: 1-3 curves → only linechart (no Multiaxis & Facet) | PASS | 0s | N/A | Confirmed: PK model (2 curves: depot, centr) shows line chart only — no Multiaxis/Facet |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 43s |
| Spec file generation | 3s |
| Spec script execution | 24s |

## Summary

All 3 functional steps passed. The pk.ivp file opened from the Files browser as a DiffStudio view. Modifying Step (0.01→0.1) and Count (1→4) updated the URL query parameters and regenerated the charts. Opening the copied URL in a new tab reproduced the exact same model state. The remark about Multiaxis/Facet absence is confirmed.

## Retrospective

### What worked well
- URL navigation to `/files/system.appdata/DiffStudio/library` reliably shows library files
- Clicking `pk.ivp` opens DiffStudio with PK model automatically
- URL parameters (`params:step=0.1&count=4`) encode all model inputs — sharing works by copying the URL
- New tab loaded same model state without friction (session preserved)
- Input modification via UI (Ctrl+A, type value, Enter) worked reliably

### What did not work
- Initial URL navigation with capitalized `Library` failed (path is case-sensitive on Linux); needed JS API `grok.dapi.files.list()` to discover correct lowercase `library`
- Double-clicking the Files tree node in Browse didn't expand child nodes; direct URL navigation was used instead

### Suggestions for the platform
- File browser URL routing could be case-insensitive for better usability
- A "Copy link" button on the DiffStudio ribbon would make sharing more discoverable

### Suggestions for the scenario
- Step 1 says "Library" (capitalized) but the actual folder is lowercase `library` — update to match
- Step 2 says "slider" and "clicker" but pk.ivp inputs appear as plain text fields — clarify control types or update model annotations
- Step 3 could note that parameters are encoded in the URL query string
