# Highlight — Run Results

**Date**: 2026-04-22
**URL**: http://localhost:8888 (local)
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog.csv dataset | 13s | PASS | PASSED | `grok.dapi.files.readCsv` + `addTableView`; 5850 rows, 11 columns. |
| 2 | Open Add New Column dialog | 26s | PASS | PASSED | Clicked `[name="icon-add-new-column"]`; dialog `[name="dialog-Add-New-Column"]` opened. Expression editor is CodeMirror (`.cm-content`). |
| 3a | Type `Abs(${AGE})` — column name highlighted blue | 1m 7s | PASS | PASSED | `${AGE}` wrapped in `span.cm-column-name`, computed color `rgb(80, 169, 197)`. |
| 3b | Type `Avg($[AGE])` — column name highlighted blue | 38s | PASS | PASSED | `$[AGE]` wrapped in `span.cm-column-name`, same blue color. |
| 4 | Add `Sqrt(` + column via `$` autocomplete — highlighted blue | 1m 37s | PASS | PASSED | MCP: typed `$`, popup shown, typed `HE`, Enter inserted `${HEIGHT}` — highlighted blue. Playwright replay initially failed (Enter closed dialog via `closeOnEnter`; and `.click()` raced CM's stability check). Fix: dispatch `mousedown` on the `li[aria-selected=true]` in the `.cm-tooltip-autocomplete` popup (CM's autocomplete listens for `mousedown`, not `click`), and accept the multi-span tokenization as well as the consolidated `${HEIGHT}` span. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~3m 57s |
| grok-browser execution (scenario steps) | ~46s |
| Execute via grok-browser (total) | 4m 43s |
| Spec file generation | 1m 13s |
| Spec script execution | 25s |
| **Total scenario run (with model)** | ~8m 20s |

## Summary

All 5 steps PASS in both the MCP run and the Playwright replay. `${AGE}`, `$[AGE]` and the autocomplete-inserted `${HEIGHT}` are all wrapped in `span.cm-column-name` with color `rgb(80, 169, 197)`. The feature under test works. Step 4 initially failed in Playwright replay — two root causes: pressing Enter in the autocomplete popup closes the dialog via `closeOnEnter`, and `.click()` on the popup `li` races CM's stability check. Replaced with a direct `mousedown` dispatch on the selected `li`; all 5 steps now PASSED.

## Retrospective

### What worked well
- `${NAME}` and `$[NAME]` syntaxes are both recognized and styled identically (`span.cm-column-name`, blue `rgb(80, 169, 197)`).
- Autocomplete-inserted column references receive the same styling as typed ones — single code path.
- The CodeMirror `.cm-column-name` class gives automation a stable semantic hook.
- Real-time re-tokenization while typing: spans update on every keystroke.

### What did not work
- Initial Playwright replay: pressing Enter after filtering the autocomplete popup to a single match closed the dialog (`closeOnEnter`) before the span could be queried. Replaced with a click on the popup `li` — that also failed because CM autocomplete listens for `mousedown`, not `click`, and Playwright's `.click()` races CM's stability check. Final fix: dispatch `mousedown` directly on the `li[aria-selected=true]` inside `.cm-tooltip-autocomplete`.
- Closing the dialog via `[name="button-CANCEL"]` was flaky at cleanup — first attempt failed silently; only `.click()` on a `uid` from a fresh snapshot worked.

### Suggestions for the platform
- When the CodeMirror autocomplete popup is open, Enter should be consumed by the popup unconditionally; it should not bubble to the dialog's `closeOnEnter`. Today, if the popup auto-dismisses after the user's last keystroke narrows the list to a single match, Enter hits the dialog instead of inserting the completion.
- Column references get retokenized into different span shapes depending on typing state (`${HEIGHT}` as one span vs. split into `$` + `{` + `HEIGHT` + `}`). Consolidating to a single span as soon as the reference is syntactically complete would make automation assertions simpler and the highlight more stable visually.
- Consider exposing a `data-column-name="HEIGHT"` attribute on `.cm-column-name` spans — cleaner for automation and would let the platform differentiate resolved vs. unresolved refs.

### Suggestions for the scenario
- Step 3 is numbered twice in `highlight.md` ("3. Paste 'Abs(...)'" and "3. Paste 'Avg(...)'") — renumber the second one to 4 and the last to 5 so the scenario has a single ordering.
- The scenario writes the column name as lowercase `age`, but `demog.csv` has uppercase `AGE`. Formula parsing is case-sensitive for the resolved-column highlight; using the exact case avoids confusion. Recommend `${AGE}` / `$[AGE]`.
- The two `${...}` and `$[...]` syntaxes aren't called out as distinct in the scenario — clarifying that both are valid column-reference notations would make the assertion intent clearer.
- Add an explicit negative case: "Paste `Abs(${nonexistent})` — column name is NOT highlighted in blue (or is highlighted differently)" to test that the editor distinguishes resolved from unresolved refs.

---
{"order": 4}
