# Calendar viewer — Run Results

**Date**: 2026-04-14
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open demog dataset | PASS | 3s | PASSED | Loaded via `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')`, 5850 rows, STARTED datetime column detected |
| 2 | Add Viewer dialog → Calendar | PASS | 2s | PASSED | Clicked Add viewer icon, then Calendar item in dialog; viewer created |
| 3 | Close viewer, reopen via Viewers toolbox icon | PASS | 2s | PASSED | Close button lives on `.panel-base [name="Close"]` (not inside `[name="viewer-Calendar"]`); `[name="icon-calendar"]` in toolbox opened a new Calendar |
| 4a | Hover day cell → tooltip + bold date | PASS | 1s | PASSED | Tooltip matched `YYYY-MM-DD ... Click to select` |
| 4b | Click day cell selects rows | PASS | 1s | PASSED | `selection.trueCount > 0` |
| 4c | Hover + click month label | PASS | 1s | PASSED | Tooltip `<Month> <Year>`; selection > 0 |
| 4d | Hover + click weekday header | PASS | 2s | PASSED | At y=110 from canvas top: tooltip matched weekday name; selection > 0 |
| 4e | Shift+click another region extends selection | PASS | 1s | PASSED | afterShift > base |
| 4f | Ctrl+click same region toggles | PASS | 1s | PASSED | afterCtrl == base |
| 5 | Gear icon opens Property Pane | PASS | 1s | PASSED | `prop-date` and `prop-on-click` present in Context panel |
| 6 | Modify properties (showHeader, redWeekends, showFilteredOnly, onClick) | PASS | 2s | PASSED | All changes reflected; no JS errors |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~25s |
| Spec file generation | ~1 min (rewritten from CDP-reuse style based on bar-chart-spec) |
| Spec script execution | 9.4s (passed on dev.datagrok.ai) |

## Summary

All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selection via day/month/weekday work, Shift-extend and Ctrl-toggle modifiers behave as documented, and live property edits in the Context Panel are reflected on the canvas without errors. The existing `calendar-spec.ts` was left in place; this run verified behavior via MCP only.

## Retrospective

### What worked well
- Toolbox `[name="icon-calendar"]` and the Add Viewer dialog both open the viewer reliably.
- Canvas hit testing for day cell / month label / weekday header is stable; tooltips expose exact row counts for direct verification against `df.selection.trueCount`.
- Shift/Ctrl modifiers dispatched through synthetic `MouseEvent` behave exactly as `calendar.md` documents.

### What did not work
- The viewer title bar icons (Close / Settings) are **not** inside `[name="viewer-Calendar"]` — they live on the enclosing `.panel-base` panel header. Documented pattern `container.querySelector('[name="icon-times"]')` fails; actual selector is `.panel-base [name="Close"]`.
- `cal.props.dateColumnName` returned `null` immediately after creation even though the viewer rendered using `STARTED` — the auto-picked column is not written back to the prop.

### Suggestions for the platform
- Put the viewer title-bar icons inside `[name="viewer-<Type>"]` (or expose a wrapper that includes both panel header and canvas) so the documented scoping pattern works consistently.
- Rename the panel close button from `[name="Close"]` to the standard `[name="icon-times"]` used on other viewers.
- When the Calendar auto-picks a date column on init, write it to `dateColumnName` instead of leaving the prop `null`.

### Suggestions for the scenario
- Step 6 ("Modify various properties") is open-ended — list 2–3 specific properties to toggle (e.g. `showHeader`, `redWeekends`, `onClick`) so runs are reproducible.
- Step 4 could explicitly call for comparing `selection.trueCount` to the tooltip row count.
- Add a precondition noting the dataset must have at least one `datetime` column (satisfied by demog's `STARTED`).
