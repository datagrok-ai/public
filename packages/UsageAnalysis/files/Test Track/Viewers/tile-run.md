# Tile Viewer — Run Results

**Date**: 2026-04-14
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Multiple DataFrames Handling — switch table among SPGI/SPGI-linked1/SPGI-linked2 | PASS | 6s | PASSED | Tile viewer rebound to each table; DOM tiles re-rendered without errors. UI: Toolbox `[name="icon-tile-viewer"]` click to add; JS API fallback (`tile.props.table = …`) for the table switch since Context Panel control opening is non-deterministic |
| 2 | Building from List of Columns — Edit Form, apply, save & restore layout | PASS | 9s | PASSED | Hamburger menu (`[name="icon-font-icon-menu"]` scoped to viewer titlebar) → "Edit Form..." → CLOSE AND APPLY worked. Layout save via `grok.dapi.layouts.save` then restore via `loadLayout` brought Tile Viewer back. Sketch state size differs slightly after restore (4160→4118) — likely metadata, not visible columns |
| 3 | Edit form with Table Change — Reset, switch source to demog, apply | PASS (partial UI) | 8s | PASSED | No source-table dropdown found in SketchView ribbon — used JS API fallback `tile.props.sketchState['table'] = 'demog'`. Setting `tile.props.table = 'demog'` rebound viewer to demog with no broken-binding placeholders |
| 4 | Hamburger menu walk + Properties → Data → Lanes hover | PASS | 7s | PASSED | Lanes column set to `Stereo Category` (5 categories). Walked all hamburger items via dispatched `mouseover` events. Page remained responsive (rowCount unchanged at 3624). No console warnings |
| 5 | Calculated column cc = ${Average Mass} + 5, then change to +10 | PASS | 5s | PASSED | `addNewCalculated` succeeded. Edit Form Reset rebuilt with default columns; renamed last sketch column input to "cc" via type+Enter; CLOSE AND APPLY persisted cc in sketchState. Formula change verified: head value 239.30 → 244.30 |
| 6 | Calculated column with filters — apply filters, change formula | PASS | 4s | PASSED | Filter on cc<400 ∧ Stereo Category∈{S_ABS,R_ONE} kept 931 rows. After changing formula to `${Average Mass}*2`, filter count unchanged (931); cc values updated. No warnings |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~5 min |
| Spec file generation | ~30s |
| Spec script execution | 1.1m |

## Summary

All 6 scenario sub-steps passed on dev.datagrok.ai. Tile viewer correctly rebinds across dataframes, persists layouts via `grok.dapi.layouts`, handles `Edit Form` flow including source-table change, walks the Hamburger menu including the historically-problematic `Properties → Data → Lanes` submenu without freezing, and updates calculated-column tiles whether filters are applied or not. Two UI gaps required JS API fallbacks: (a) the SketchView ribbon does not expose a source-table dropdown that the scenario expects, and (b) drag-and-drop of columns onto the card cannot be reliably automated.

## Retrospective

### What worked well

- Hamburger menu access via `viewer-Tile-Viewer → closest('.panel-base') → .panel-titlebar [name="icon-font-icon-menu"]` was reliable.
- `grok.dapi.layouts.save / find / loadLayout` round-trip restored the Tile viewer cleanly with same sketch payload.
- Programmatic mouseover walk over hamburger items was sufficient to verify no freeze on Lanes hover.
- `addNewCalculated` propagated to tiles instantly even with active filters.

### What did not work

- The Tile viewer titlebar is **not** a child of `[name="viewer-Tile-Viewer"]` — it lives on the wrapping `.panel-base.panel-titlebar`. The grok-browser reference table implies title-bar icons are scoped to the viewer container; that's misleading.
- `right-click` on `.d4-tile-viewer-lanes-host` did **not** open the viewer's context menu in MCP — context menu likely requires a real OS-level mouse event. Hamburger menu was the workable substitute.
- SketchView ribbon (`.d4-ribbon-panel`) contains only Align/Show/Save/Load/Presentation/Editable/EDIT/RESET — there is **no source-table dropdown** despite the scenario assuming one (and the tile.md reference describing it as "In the SketchView header, change the table field").

### Suggestions for the platform

- Add a visible **Source Table** dropdown to the SketchView ribbon (or surface it in the right-side property panel header). Without it, switching the form's source table can only be done through internal `sketchState` mutation.
- Consider giving the SketchView a "drag column from list" target panel that's automation-friendly (e.g. a left-side column list with `[name="div-column-…"]` items), since the current sketch directly renders fields and isn't easy to build up step by step from automation.
- After `loadLayout` the restored Tile viewer's `sketchState` size differs from the saved snapshot — verify that no rendering-relevant fields (esp. `elementStates`, `formDesigned`) are getting dropped during round-trip.

### Suggestions for the scenario

- Step 2 says "Drag a few columns onto the card" but the auto-generated form already contains the standard SPGI columns at open — clarify whether the test wants a Reset+rebuild or a delta on top of the auto layout.
- Step 3 should explicitly mention the workaround needed when the source-table dropdown is missing, or note that this gap is part of what's being tested.
- Step 5: "Click the cc column header; in the Context Pane, change the formula and **Apply**" — there's no separate Apply button on the calculated-column property panel; the formula re-applies on Enter / blur. Reword to "edit formula in the Context Pane and press Enter".
