# Formula Refreshing — Run Results

**Date**: 2026-04-22
**URL**: http://localhost:8888
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog.csv | 4s | PASS | PASSED | 5850 rows, 11 columns, WEIGHT present. Loaded via `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')`. |
| 2 | Create Weight2 = `${WEIGHT} + 100` | 35s | PASS | PASSED | Toolbar `[name="icon-add-new-column"]` opens dialog. Row 0: WEIGHT=73.2, Weight2=173.2. Formula tag recorded. |
| 3 | Create Weight3 = `${Weight2} + 100` | 15s | PASS | PASSED | Row 0: Weight3=273.2 = Weight2+100. |
| 4 | Create Weight4 = `Log10(${Weight3}) - 0.2` | 15s | PASS | PASSED | Row 0: Weight4=2.2365. All 5 sampled rows matched. |
| 5a | Edit Weight2 formula → propagate to Weight3, Weight4 | 20s | PASS | PASSED | Formula pane `.d4-pane-formula` used as selector. `${WEIGHT} + 200` applied; Weight2=273.2, Weight3=373.2, Weight4=2.372. |
| 5b | Edit Weight3 formula → propagate to Weight4 | 15s | PASS | PASSED | `${Weight2} * 2` applied → Weight3=546.4, Weight4=2.538. |
| 5c | Edit Weight4 formula | 15s | PASS | PASSED | `Log10(${Weight3}) + 1` applied → Weight4=3.738. |

**Time** = wall-clock per step. **Result** = MCP outcome. **Playwright** = spec outcome (after fixes).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~2m |
| grok-browser execution (scenario steps) | ~15s |
| Execute via grok-browser (total) | ~2m 15s |
| Spec file generation | ~1m |
| Spec script execution | 28s |
| **Total scenario run (with model)** | ~3m 45s |

## Summary

All seven sub-steps pass in both the MCP run and the Playwright replay. Dependency propagation across calculated columns is verified end-to-end: Weight2 → Weight3 → Weight4 chain recomputes correctly after in-place formula edits via the Context Panel's Formula accordion. Two repeat-each runs pass cleanly.

## Playwright fixes applied

The prior Playwright replay failed with three distinct issues; all are now resolved:

1. **Cross-talk with legacy spec** — `/home/andrew/tmp-pw/anc-formula-refreshing.spec.ts` (a legacy variant covering the WEIGHT[0]=200 re-fire scenario) ran in parallel against the same Datagrok server and corrupted shared state. Archived to `anc-formula-refreshing.spec.ts.legacy` so Playwright no longer picks it up.
2. **Preloader overlay intercepting clicks** — `<div id="grok-preloader" class="grok-preloader">` (z-index 100500) occasionally covers the dialog OK button right after a calc-column add. Playwright's `click()` (even with `force: true`) respects browser hit-testing, so it kept retrying. Fix: dispatch MouseEvents (`mousedown`/`mouseup`/`click`) directly on the target via `page.evaluate`, which bypasses hit-testing entirely.
3. **Stale Formula-pane CodeMirror** — after `g.shell.o = col`, the `.d4-pane-formula .cm-content` takes ~500–1500 ms to refresh with the newly-selected column's formula. The prior spec only waited for the pane element to exist, so it would type `${Weight2} * 2` into the previous column's CM (showing `${WEIGHT} + 200`), and Apply would silently keep the original formula. Fix: read the target column's current `formula` tag first, then `waitForFunction` until `.d4-pane-formula .cm-content` matches it — this confirms the pane has rebound to the new column before typing.

Secondary hardening: used `.d4-pane-formula` (stable CSS class on the pane body) instead of iterating `.d4-accordion-pane` and matching header text "Formula"; scoped the Add-New-Column `[name="button-Apply"]` click to `.d4-pane-formula [name="button-Apply"]` so it can't collide with other Apply buttons; scoped dialog lookups to `.last()` since a stuck dialog from a prior step can produce duplicates.

## Retrospective

### What worked well
- Dependency propagation through calculated-column tags is correct for Weight2→Weight3→Weight4 chain (all rows match within 1e-3).
- The Formula accordion's `.d4-pane-formula` class is stable and discoverable — much more reliable than matching accordion header text.
- `col.getTag('formula')` exposes the formula string for both read (pre-select check) and assertion.

### What did not work
- **Accordion "expanded" class lags the actual content state** — `.d4-accordion-pane-expanded` is set eventually but not synchronously with the pane rebuild. Relying on class toggling for readiness is unsafe; content-match is the correct signal.
- **`page.click({force: true})` does not bypass hit-testing overlays** — `force: true` only skips actionability checks (visible/enabled/stable). Playwright still sends a real pointer event through the browser's hit-testing, which routes the click to whatever element is on top (the preloader overlay, here). Direct `dispatchEvent` is the correct workaround.

### Suggestions for the platform
- Provide a stable public JS API to mutate a calculated column's formula and trigger dependent-column recomputation, e.g. `col.setFormula(expr)` or `df.columns.updateFormula(name, expr)`. The current only-programmable path is through the Formula pane's Apply button or the Edit-in-dialog flow — neither is easy to drive from automation.
- Emit a stable `formula-applied` event on `grok.shell.tv.dataFrame` (or on the column) so tests can `await` instead of `waitForTimeout(1500)`.
- Surface pane-body-rebuild completion via an authoritative signal (e.g. a `data-formula-for="Weight3"` attribute on the CodeMirror, or an event on `grok.shell.onCurrentObjectChanged` that fires after the context panel finishes rendering) — would eliminate the content-match heuristic.
- Add stable, globally-unique `name` attributes to per-pane Apply/Edit-in-dialog buttons (e.g. `button-formula-pane-Apply`) — the bare `button-Apply` selector is ambiguous across panels.
- Consider scoping the `#grok-preloader` element or its z-index so it never covers open dialogs — the current overlay behavior caused several minutes of debugging.

### Suggestions for the scenario
- Step 5 currently says "Open the Context panel. Modify the formulas for Weight2, Weight3, and Weight4" — spell out the exact path (select column → Context Panel → Formula accordion → edit CM → Apply) so automation doesn't have to guess between Context-Panel-Apply vs. right-click-Edit-Formula-dialog.
- Add explicit numeric post-conditions after step 5 (e.g. "After setting Weight2 formula to `${WEIGHT} + 200`, Weight3 row 0 should be `WEIGHT[0] + 300`, Weight4 row 0 should be `Log10(WEIGHT[0] + 300) - 0.2`"). Exact expected values let automated checks avoid subjective "looks right" judgement.
- Split "Additional Notes: save as project and re-open" into its own scenario — project serialization + reopen exercises a separate code path from in-memory dependency recalc.


---
{"order": 7}
