# Calendar (Playwright) — Run Results

**Date**: 2026-08-06
**URL**: https://dev.datagrok.ai (1.28.0)
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Calendar from the Viewers toolbox | PASS | Canvas carries content; `STARTED` picked up automatically |
| 2 | Day, month and weekday tooltips report the matching rows | PASS | Day `1989-12-06` → 6 rows, `December 1989` → 264 rows, weekdays Mon–Sat match the data; Sunday is labelled correctly but reports 0 — `[KNOWN_BUG_REPRODUCED:GROK-20634]` |
| 3 | The auto-picked date column is reflected in the property grid | PASS | `dateColumnName` is `null` while `STARTED` is drawn — `[KNOWN_BUG_REPRODUCED:GROK-20636]` |
| 4 | Clicking a day, a month and a weekday selects their rows | PASS | Selection equals the tooltip's count in all three regions |
| 5 | Shift+click extends the selection, Ctrl+click toggles it | PASS | Monday 849 → +Tuesday 837 → back to 849 |
| 6 | On Click set to Filter makes the same click filter instead | PASS | Tooltip switches to *Click to filter*; a month click filters to 264 rows and selects nothing |
| 7 | Red Weekends colours the weekend day numbers | PASS | 2273 red pixels → 0 → back |
| 8 | Show Header hides the year caption | PASS | Content pixels 345832 → 392045; the weekday row hit-tests correctly at its new position |
| 9 | Filtering the table reshapes the calendar | PASS | Canvas repaints; tooltips keep counting all rows behind a date; toggling **Show Filtered Only** under the filter changes nothing — `[KNOWN_BUG_REPRODUCED:GROK-20635]` |
| 10 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 30.3s / 30.5s on two consecutive runs (was 19.3s before the three guarded steps, which spend their timeouts waiting for a change that never comes) |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `expect.poll`).

* The calendar exposes nothing inside its canvas, so the three hit-test regions are
  derived from the same geometry the viewer uses: `dy = (width - 5) / 8`, a month
  column `dy + 5` wide on the left, a day-of-week row `dy` high under an optional
  year header of `dy * 1.1`, and the day cells below. Passing `withHeader: false`
  to the same helper is what lets step 7 prove the row really moved up.
* Day cells are not addressed by date — the spec probes a few cells and works with
  the first one that carries rows, then checks its count against the data. That
  keeps the step independent of where demog's dates happen to fall.
* Red weekends are asserted with a red `rgbRange` on `countCanvasPixels`, which is
  exact (2273 → 0) rather than a "something repainted" delta.

## Open bugs — guarded by `knownOpenBug()`

| Bug | Detail |
|---|---|
| GROK-20634 — Sunday column selects nothing | The hit test compares `date.weekday == dayOfWeek`; the column index is 0-based with Sunday first, Dart's `weekday` is 1..7 with Sunday = 7. Sunday reports 0 rows and selects nothing — 837 rows on demog, confirmed by counting the `STARTED` column directly. Monday–Saturday are asserted hard, as is Sunday's label; only its count is guarded |
| GROK-20635 — Show Filtered Only does nothing | `showFilteredOnly` is read only by `CalendarCore.rowIndexes`, which is never called; `_refresh()` walks every row regardless. Asserted under an active `SEX = M` filter, where the setting would have to matter |
| GROK-20636 — Date selector not wired up | The viewer renders from `STARTED` while `dateColumnName` stays `null`, and nothing reads that property back either, so picking a column does not change what is drawn. demog has one date column, so only the write-back half is asserted here; the other half needs `System:DemoFiles/chem/SPGI.csv`, which has several date columns with very different ranges |

## Not automated

* Adding the viewer from the **Add viewer** dialog, closing and reopening it from
  the toolbox, the bold-on-hover day number, and the *Style* colours — kept as a
  manual list in `calendar.md`.
