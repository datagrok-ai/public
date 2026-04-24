# Composition Analysis manual test — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open sample_FASTA / sample_HELM / sample_MSA | 45s | PASS | PASSED | Opened `System:AppData/Bio/samples/{FASTA,HELM,MSA}.csv` via `grok.dapi.files.readCsv`; Macromolecule semType detected on Sequence/HELM/MSA columns (64/540/540 rows). |
| 2 | Bio > Analyze > Composition opens WebLogo viewer | 30s | PASS | PASSED | Menu dispatched via `page.evaluate` (click `[name="div-Bio"]` → `mouseenter` on `[name="div-Bio---Analyze"]` → click `[name="div-Bio---Analyze---Composition"]`). `page.locator(...).click()` fails because the top-menu item is in DOM but not visible until a parent submenu is opened; the evaluate-based dispatch bypasses Playwright's visibility gate. Viewer presence verified via `waitForFunction` on `grok.shell.tv.viewers`, not a DOM locator. Two toast errors still appear on FASTA first-attach: `Input not found for property type "object" for property "rules"` and `TypeError: Cannot read properties of null (reading 'length')`. |
| 3 | Click a letter in the WebLogo — rows get selected | 20s | PASS | PASSED | After WebLogo is added to the viewers list we still need a ~3s settle before the canvas is interactive — without it, the click fires but no rows get selected on the 2nd/3rd dataset. Click dispatched as `mousedown`/`mouseup`/`click` `MouseEvent` against the canvas at multiple x-offsets (`18, 30, 50, 80, 120, 160`) so narrower columns on HELM/MSA still land on a letter cell. First offset to select ≥1 row wins. Selection counts: FASTA 47/64, HELM 435/540, MSA 435/540. |
| 4 | Click the Gear icon on top of the viewer | 25s | PASS | PASSED | Gear lives in the outer `.panel-titlebar [name="icon-font-icon-settings"]`, NOT inside the `[name="viewer-WebLogo"]` subtree — scope via the viewer's `.panel-base` ancestor. Gear click opens `.grok-prop-panel .property-grid`. `tr[name="prop-show-position-labels"]` exists but is hidden (inside a collapsed-by-default category), so wait with `state: 'attached'` rather than the default `visible`. |
| 5 | Go to Context Pane and change arbitrary properties | 25s | PASS | PASSED | Toggled `Show Position Labels` and `Skip Empty Positions` checkboxes via `tr[name="prop-show-position-labels"]` / `tr[name="prop-skip-empty-positions"]`. At least one of the two boolean properties flipped its value for each dataset (confirmed by diffing `wl.getOptions().look.*` before/after). |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m 5s |
| grok-browser execution (scenario steps) | 2m 8s |
| Execute via grok-browser (total) | 4m 13s |
| Spec file generation | 1m 29s |
| Spec script execution | 41s |
| **Total scenario run (with model)** | 6m 23s |

## Summary

All 15 checks pass in both browser and Playwright (5 scenario steps × 3 datasets): WebLogo opens from **Bio → Analyze → Composition**, clicking a letter selects rows (47/64 on FASTA, 435/540 on HELM and MSA), the gear icon opens the property grid, and at least one property flips via the Context Pane for each dataset. Spec run dropped from 1m 6s to 41s after three fixes (see Retrospective). Total scenario run (with model): 6m 23s.

## Retrospective

### What worked well
- `grok.dapi.files.readCsv` + `addTableView` + `onSemanticTypeDetected` + a 5s Bio/Chem settle reliably brings up a Bio dataset with Macromolecule rendering.
- The top-menu item for `Bio > Analyze > Composition` has a stable `[name="div-Bio---Analyze---Composition"]` attribute with no `...` suffix (no dialog). Dispatching `.click()` inside `page.evaluate` invokes the menu action directly and works identically in Playwright.
- A simple `MouseEvent` click at canvas coordinates drives WebLogo letter selection; probing a short list of x-offsets makes this robust across FASTA/HELM/MSA where the letter-column pitch differs.
- Scoping `.panel-base > .panel-titlebar [name="icon-font-icon-settings"]` is the cleanest way to find a docked viewer's gear when the viewer's own root doesn't host a title bar.
- Property-grid rows use the stable `tr[name="prop-<kebab>"]` naming; boolean properties expose an inner `input[type="checkbox"]` that responds to `cb.click()` + `change` event.

### What did not work initially (and how it was fixed)
- `page.locator('[name="div-Bio---Analyze---Composition"]').click()` — Playwright enforces visibility and the menu item is in DOM but hidden. Fixed by dispatching the click via `page.evaluate` (click parent `div-Bio` → `mouseenter` on `Analyze` → click `Composition`).
- `page.locator('[name="viewer-WebLogo"]').waitFor(...)` — the `[name="viewer-WebLogo"]` wrapper lives on the panel root but may not be `visible` under Playwright's rules during a docking transition. Fixed by waiting on `grok.shell.tv.viewers` via `waitForFunction`.
- Canvas click returned 0 selected rows on the 2nd/3rd dataset — the canvas was in DOM but WebLogo hadn't fully rendered/wired handlers yet. Fixed with a 3s settle after step 2 plus probing multiple x-offsets (`18, 30, 50, 80, 120, 160`).
- `tr[name="prop-show-position-labels"]` was hidden (collapsed accordion section) so the default `visible` wait timed out. Fixed with `state: 'attached'`.

### What did not work (platform issues, not fixed)
- The right-side property sidebar shows two distracting toast errors the first time a WebLogo is attached on FASTA: `Input not found for property type "object" for property "rules"` and `TypeError: Cannot read properties of null (reading 'length')`. They don't block interaction but are visible to users.
- Numeric property cells (`Position Width`, `Max Monomer Letters`) did not commit new values with a simple click → select → `insertText` + Enter sequence; only boolean checkbox toggles flipped state cleanly. The editable `<input>` isn't `document.activeElement` after the cell click.

### Suggestions for the platform
- Put `name="icon-font-icon-settings"` (and the rest of the title-bar icon set) **inside** the WebLogo viewer's own root so that test selectors scoped to `[name="viewer-WebLogo"]` find the gear without walking up to the docking panel.
- Fix the `rules` property-type resolution and the `Cannot read properties of null (reading 'length')` error that appears in the right-side panel when a WebLogo is first attached.
- Make numeric property-grid cells enter edit mode on a single click (focus the `<input>` reliably) so keyboard-driven automation matches user expectation.
- Surface top-menu items in a way that passes Playwright's visibility gate from a cold page — e.g. an API like `grok.shell.topMenu.find('Bio/Analyze/Composition').click()`, or an ARIA-exposed toggle that opens the submenu first. The evaluate-and-click workaround replicates across every spec that drives the top menu.
- Expand property-grid accordion sections by default (or expose a stable `aria-expanded` attribute) so automation can assert property rows without special-casing hidden vs. visible waits.
- Emit a viewer-ready event after render+bind is complete (beyond `viewers.find(...).type === 'WebLogo'` showing up in the list). Currently every spec has to add ~3s of wall-clock guessing for canvas interactivity.

### Suggestions for the scenario
- Replace the filenames `sample_FASTA.csv`, `sample_HELM.csv`, `sample_MSA.csv` with the actual paths: `System:AppData/Bio/samples/FASTA.csv`, `.../HELM.csv`, `.../MSA.csv` — no file named `sample_*` exists in that folder.
- Spell out what "change arbitrary properties" means in terms of expected outcome — e.g. "toggle Show Position Labels, confirm labels disappear under the logo" — so the tester has a clear PASS/FAIL.
- Add a final step: close the WebLogo viewer (right-click → Close or `Ctrl+F4`) and confirm the grid view is restored — catches the common dirty-state issue when chaining scenarios.
