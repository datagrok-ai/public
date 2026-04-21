# PLS — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open cars.csv from Demo files | 1s | PASS | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/cars.csv')` + `addTableView`; 30 rows, 17 cols; `onSemanticTypeDetected` resolved in ~1s. Selenium/tabs/showFilters flags set in one batch |
| 2 | Top Menu > ML > Analyze > PLS... | 2s | PASS | PASSED | `[name="div-ML---Analyze---PLS..."]` clicked via `document.querySelector(...).click()`; `.d4-dialog` opened with title "PLS"; inputs `input-host-Predict/Using/Components/Quadratic` present (Names hidden). Defaults: Predict=price, Using empty, Components=2 |
| 3 | Select all Using features; set Components = 3 | 4s | PARTIAL | PASSED | `[name="input-host-Using"] .ui-input-editor` opened "Select columns..."; `[name="label-All"]` + top-most `[name="button-OK"]` selected 16 numerics (string `model` excluded); Using display = "(16) All". `[name="input-Components"]` + Control+A + type "3" -> DOM value "3". Validation conflict surfaced: `input-host-Predict` turns red because `price` is also in Using, and RUN stays disabled |
| 4 | Execute PLS via RUN; expect PLS1/PLS2/PLS3 | 1m 30s | FAIL | FAILED | `[name="button-RUN"]` stays `disabled`; clicking does nothing; no columns added. Spec asserted disabled-state explicitly and threw with a diagnostic message (Predict + Using hosts have `color: rgb(235, 103, 103); border-bottom-color: rgb(235, 103, 103)`). Scenario wording "select all available columns" conflicts with PLS validation (Predict must be excluded from Using); no tooltip/help on disabled RUN explains why |
| 4b (spec-only) | Fallback — deselect `price` from Using and retry | n/a | n/a | FAILED | Spec reopened "Select columns..." and tried to toggle `price` off by clicking a label with that text; toggle did not register. `[name="label-price"]` fallback was not present. Kept as diagnostic path to prove the underlying issue is input-validation UX, not WASM/server |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 30s |
| grok-browser execution (scenario steps) | 10s |
| Execute via grok-browser (total) | 2m |
| Spec file generation | 2m |
| Spec script execution | 15s |
| **Total scenario run (with model)** | 4m 15s |

## Summary

MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 2 PASS / 1 PARTIAL / 1 FAIL. The
dialog path (menu, Using selector, Components) works, but PLS validation disables RUN when the
Predict column is also in Using — and the scenario's "select all available columns" wording leaves
`price` in both lists. Nothing in the UI tells the user why RUN is disabled; the only cue is a red
border on Predict. The generated Playwright spec reproduced the disabled-RUN failure cleanly
(step 4 threw with diagnostic styles), and the fallback step that attempts to deselect `price`
from Using also failed — the toggle path via label click didn't register, which is itself a minor
automation gap in the "Select columns..." sub-dialog. **Total scenario run (with model): 4m 15s.**

## Retrospective

### What worked well
- `[name="div-ML---Analyze---PLS..."]` is stable for the PLS menu entry; `document.querySelector(...).click()` is enough — no hover-chain needed.
- `[name="input-host-Using"] .ui-input-editor` + `[name="label-All"]` + top-most `[name="button-OK"]` cleanly picks all 16 numeric columns.
- `[name="input-Components"]` accepts real click + `Control+A` + type "3" + `Tab` and the DOM value updates to "3".
- Polling `grok.shell.tv.dataFrame.columns` names is a deterministic "did PLS produce PLS1/2/3" signal.
- Explicit RUN-enabled check in the spec caught the disabled state with a precise diagnostic (input border colors), making the platform issue unambiguous rather than a generic "click did nothing".

### What did not work
- PLS validation disables RUN when the Predict column appears in Using, but there is no tooltip, banner, or status line on the disabled RUN button explaining why. The only indicator is a red border on the Predict selector, which is silent and easy to miss.
- The "select all available columns" scenario wording directly contradicts PLS validation: after clicking All, `price` is in both Predict and Using, so RUN is immediately disabled. A user following the scenario literally cannot proceed.
- The fallback toggle-off-`price` step in the spec did not register: clicking an element whose `textContent` is `price` inside the "Select columns..." sub-dialog did not toggle the checkbox, and `[name="label-price"]` was not present. The column list in this sub-dialog does not expose per-column `name` attributes, making deterministic automation of individual rows hard.

### Suggestions for the platform
- **When RUN is disabled due to input validation, surface the reason inline (tooltip or banner) so the user understands why. Red border on Predict alone isn't discoverable.** Preferred: hover-tooltip on the disabled RUN button ("RUN disabled: Predict column 'price' must be removed from Using") and/or a small banner under the dialog title.
- Auto-exclude the Predict column from the Using multi-select (either by filtering its options, or by automatically deselecting it when the user chooses All) so "select all" does the expected thing.
- Add per-column `name` attributes (`[name="label-<col>"]` or `[name="row-<col>"]`) to the "Select columns..." sub-dialog rows so automation and assistive tech can address individual entries.
- Consider showing an inline error message under the Predict field ("Cannot be in Using") rather than only a color change — color-alone fails WCAG 1.4.1 (use of color).

### Suggestions for the scenario
- Clarify step 3: "Select all available columns **except the Predict column**" — current wording causes a validation conflict unrelated to the PLS behaviour under test.
- Add an explicit sub-assertion for step 4: "verify `df.columns.names()` now contains PLS1, PLS2, PLS3" so a tester catches silent failure modes.
- Add a secondary step that explicitly tests the validation UX: "set Using to include Predict; verify RUN is disabled and an error/tooltip explains why".

### Suggestions for the skill
- Document that Playwright's default `testMatch` rejects `-spec.ts`; either adopt `.spec.ts` naming or document the `--config /tmp/pw-config-*.ts` workaround alongside the run instructions (recurring issue across multiple prior runs).
- When a scenario is expected to expose a validation-UX issue, recommend splitting the spec into (a) the disabled-RUN assertion and (b) a fallback that proves the tool works with valid input — the pattern used here.
