# Add and Edit Sticky Metadata — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | S1 Open SPGI dataset | 10s | PASS | PASSED | 3624 rows, 88 cols (fresh); Structure semType=Molecule |
| 2 | S1 Right-click Structure cell r0 → Sticky meta → Edit for current cell | 8s | PASS | PASSED | Dialog "Sticky meta" opened with TestSchema1 section (rating/notes/verified/review_date/approve) |
| 3 | S1 Fill rating=5, notes="test note", verified=true, save | 15s | PASS | PASSED | review_date auto-populated (4/23/2026 7:38 AM); first button-Save clicked for TestSchema1 |
| 4 | S1 Close + re-open dialog, verify metadata persisted | 8s | PASS | PASSED | On re-open: rating=5, notes="test note", verified=true retained |
| 5 | S2 Add TestSchema1 sticky columns via Context Panel (+ all) | 6s | PASS | PASSED | 5 sticky cols added (rating, notes, verified, review_date, approve). Schema-match add is async — poll up to 15s for cols to appear |
| 6 | S2 Sort ascending by rating column | 3s | PASS | PASSED | grid.sort(['rating'], [true]); ascending arrow visible; 3 → 5 → 55… |
| 7 | S2 Remove rating column and re-add; metadata persists | 6s | PASS | PASSED | Removal via df.columns.remove; re-add via + button in TestSchema1 section; rating=5 preserved for row CAST=634783 |
| 8 | S3 Select 3 rows, right-click → Edit for all properties | 7s | PASS | PASSED | page.mouse right-click on Structure cell → DOM-click "Edit for all properties" (hidden menu item); batch-edit host appears |
| 9 | S3 Set verified=true and notes="batch note" for selected rows | 8s | PASS | PASSED | prop-view-verified checkbox toggled via locator; notes edit-label click → Ctrl+A/Delete/type "batch note"; all 3 rows updated |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~4m |
| grok-browser execution (scenario steps) | ~1m 11s |
| Execute via grok-browser (total) | ~5m 11s |
| Spec file generation | ~1m |
| Spec script execution | 44s |
| **Total scenario run (with model)** | ~7m |

## Summary

All 3 scenario sections (single-cell add/edit, sticky column behavior, batch edit) were reproduced successfully via MCP browser on dev.datagrok.ai as user `agolovko` and subsequently replayed end-to-end in Playwright (9 of 9 steps PASSED, 44s wall-clock). Metadata values persisted across dialog close/re-open and across column remove/re-add. Total scenario run ≈ 7 minutes.

## Retrospective

### What worked well
- Right-click → contextmenu dispatch with DOM `MouseEvent('contextmenu', {clientX, clientY, button:2})` reliably opens the grid cell menu.
- `input-host-Rating`, `input-host-Notes`, `input-host-Verified` selectors plus keyboard `type_text` correctly set TestSchema1 properties.
- Clicking the first `[name="button-Save"]` inside the multi-section Sticky meta dialog commits only TestSchema1; the button-color flash (rgb(80,169,197)) indicates the save fired.
- The property grid inside the "Edit for all properties" batch-edit view enters edit mode when the `.property-grid-item-view-label` inside the value cell is clicked (not the cell itself, not via double-click).
- Metadata fully persists across column remove + re-add — confirmed via CAST Idea ID `634783` re-lookup after sort.

### What did not work
- Setting `grok.shell.o = cell` did NOT populate the Sticky meta pane with editable inputs. The pane only renders "+ add as column" buttons for Column objects, not edit widgets for a cell. The real edit path is **grid cell right-click → Sticky meta → Edit for current cell…** (dialog), not the Context Panel pane — the scenario's wording misled the first attempt.
- Setting `grok.shell.o = row` yields only a "Values" pane, no Sticky meta widget.
- `button-Add-TestSchema1's-properties-as-columns` takes ~5 s of async schema matching before sticky columns appear. Immediate column-count checks after the click see 0 new columns. Fix: poll `df.columns.contains(...)` up to 15 s.
- `button-Add-<field>-as-a-column` names collide across multiple schemas (TestSchema1 and scratchSchema-apz8 both expose `button-Add-notes-as-a-column`). A global `querySelector('[name="..."]')` finds the wrong schema's button. Fix: scope queries to the schema's `.d4-build-root.ui-form` section whose `.d4-flex-row` header text equals `TestSchema1`.
- Menu item "Edit for all properties" is in DOM but visually hidden until its parent Sticky meta submenu is hovered. `locator.click()` with Playwright's visibility guard times out ("resolved to hidden"). Fix: DOM `.click()` via `page.evaluate` bypasses the visibility check and still fires the handler.
- Canvas grid right-click: `elementFromPoint(...).dispatchEvent(new MouseEvent('contextmenu', ...))` is flaky on the grid canvas (sometimes no menu appears). `page.mouse.move(x,y)` + `page.mouse.click(x, y, {button:'right'})` is reliable because it dispatches the full mousedown/mouseup/contextmenu sequence.
- Batch edit note field retained "batch note" from a prior run, so a naive `keyboard.type('batch note')` produced "batch notebatch note". Fix: `Ctrl+A` + `Delete` before typing.
- `df.rows.clearSelection()` does not exist in this JS API version — use `df.selection.setAll(false)`.
- `grok.shell.t` becomes `null` while the Edit-for-all-properties view is active; use `grok.shell.tables[0]` for read-back verification.

### Suggestions for the platform
- The Sticky meta Context Panel pane, when the current object is a **Cell** (not a Column), should show editable inputs for that cell's metadata — today it does not render anything useful for Cells. That would align the pane behavior with the scenario wording.
- The `+` buttons in the Sticky meta pane should toggle (show `−` when column already exists) so users can tell whether clicking adds or is a no-op.
- Context-menu "Edit for all properties" should remain stable against `grok.shell.t` access — currently the current-view switch nulls `grok.shell.t` while the batch-edit view is active.
- Tooltip on Sticky meta pane + icons to disambiguate "add column" vs "edit value".

### Suggestions for the scenario
- Step 1.2 ("click on cell under column Structure") and 1.3 ("Open Context Panel → Sticky Meta") are misleading — the Context Panel → Sticky meta pane does NOT show edit fields for a cell. Reword to: "Right-click the Structure cell → Sticky meta → Edit for current cell…". If Context Panel editing is the intended flow, that pane should be fixed (see platform suggestion above).
- Step 3 typo: `Edot` → `Edit`. Menu item is literally "Edit for all properties" (no "Sticky Meta >" prefix — it appears directly under the Sticky meta submenu).
- Clarify that "blue marker in top-right of cell" is a small colored dot (visible in canvas grid) — not a tooltip.
- Section 3 step 2 says "In popup menu choose Sticky Meta > Edot for all properties (Rows = Selected)". The Rows combo-box defaults to `Selected` when multiple rows are selected — make this explicit.
- Pre-condition should state that `create-schema-and-type.md` must be run first and that TestSchema1 must already have a matching expression covering the `Structure` Molecule column.
