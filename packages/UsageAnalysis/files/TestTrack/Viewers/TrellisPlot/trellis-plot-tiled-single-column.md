---
feature: trellisplot
realizes_atlas:
  - trellisplot.cp.tiled-single-column
  - trellisplot.int.tiled-single-column-layout
realizes:
  - viewers.trellis-plot
priority: p2
target_layer: playwright
coverage_type: edge
related_bugs: []
realized_as:
  - trellis-plot-tiled-single-column-spec.ts
expected_results:
  - anchor: "Scenario 1 — Step 4"
    expectation: >-
      Switching Use Tiled View back on — after the entry state was switched off,
      Tiled View being on by default — rebuilds the layout: the single untiled
      band of five cells becomes a two-row grid four cells wide, eight
      .d4-trellis-plot-cell slots in all. Six of them carry a category and are
      drawn; the two trailing slots of the last row are empty grid padding. Both
      the row count and the slot count differ from the untiled band read a
      moment earlier — that difference is the layout rebuild the scenario exists
      for.
  - anchor: "Scenario 1 — Step 5"
    expectation: >-
      At the default Tiles Per Row the widest row holds more than one cell;
      setting Tiles Per Row to 1 re-flows the cells into a single-column strip
      (one cell per row). The strip runs into the platform's five-cell window
      before it runs out of categories, so the number of rendered cells drops
      from the default grid's eight slots to five.
  - anchor: "Scenario 1 — Step 6"
    expectation: >-
      After changing Tiles Per Row to 3 the widest row holds exactly three cells
      — the split column has more categories than that, so the requested width
      is the one actually laid out: narrower than the default four, wider than
      the one-cell strip of the previous step. The grid is two rows, six cells.
  - anchor: "Scenario 1 — Step 7"
    expectation: >-
      Switching Use Tiled View off lays the cells out as a single horizontal
      band — exactly one row, every rendered cell side by side (geometry
      round-trip): five cells, with fewer rows and a different cell count than
      the two-row tiled grid read just before the toggle. The band is capped at
      the five cells the platform renders across, so the sixth category of the
      split column is not on screen.
  - anchor: "Scenario 2 — Step 3"
    expectation: >-
      With Auto Layout on and the viewer box large (height > 300, width > 100)
      the control panel node is present and displayed — it exists in the DOM and
      its computed display is not none.
  - anchor: "Scenario 2 — Step 4"
    expectation: >-
      After the viewer is resized below the threshold (height <= 300 or width <=
      100) the control panel's computed display goes to none — the platform
      keeps the node, it only stops displaying it — while look.showControlPanel
      still reads its last explicitly-set value (unchanged).
  - anchor: "Scenario 2 — Step 5"
    expectation: >-
      Restoring the viewer to a large size brings the control panel's computed
      display back out of none — the same node, still present throughout, is
      shown again — without changing look.showControlPanel.
  - anchor: "Scenario 2 — Step 6"
    expectation: >-
      With Auto Layout off, the control panel stays displayed at any viewer box
      size while showControlPanel is explicitly true; forcing showControlPanel
      to false collapses its computed display to none at any size — and here too
      the node itself is kept, only its display changes.
---

# Trellis Plot — Tiled Single-Column Layout and Auto Layout Coupling

## Setup

1. Close all open tables and viewers.
2. Open the demog golden dataset from System:DemoFiles/demog.csv using the
   Files browser (Browse > System > DemoFiles, double-click demog.csv).
3. Add a Trellis Plot to the table view: click the viewer icon in the
   toolbar and select Trellis Plot from the Data Separation group.
4. In the Properties panel (gear icon or Context Panel), set X Column to
   DIS_POP and leave Y Column empty so only one axis is split. DIS_POP has
   six categories — more than the default Tiles Per Row of 4 — which is what
   makes the tiled and untiled layouts differ visibly.

## Scenarios

### Scenario 1: Tiled single-column geometry and Tiles Per Row round-trip

Steps:

1. With one split column (DIS_POP, six categories), confirm the entry state:
   Use Tiled View is on — that is the shipped default for a single split
   column, so the ladder has to start by leaving it.
2. Open the Trellis Plot's Properties panel from the gear icon on the
   viewer title bar, so the Tiled View and Tiles Per Row controls are in
   reach.
3. Switch Use Tiled View off — the non-default side of the toggle — and
   confirm the viewer renders a single horizontal band of five cells: the
   layout shows at most five cells at a time, so the sixth category of the
   split column stays off screen.
4. Switch Use Tiled View back on; verify the band is rebuilt into a grid of
   two rows, four cells wide — eight cells in all, six of them drawn for the
   six categories and the two trailing ones left as empty grid padding. Both
   the number of rows and the number of cells differ from the untiled band,
   which is the layout change the toggle is for.
5. With Tiles Per Row left at its default, confirm the widest row holds
   more than one cell; then set Tiles Per Row to 1 and verify that all
   cells reflow into a single-column strip, one cell per row and five rows
   tall — the strip reaches the five-cell window before it runs out of
   categories, so five cells are shown instead of the default grid's eight.
6. Set Tiles Per Row to 3; verify the widest row now holds exactly three
   cells — narrower than the default four, wider than the one-cell strip of
   the previous step — and the grid is two rows of three, six cells in all.
7. Restore Tiles Per Row to its default, then set Use Tiled View to off;
   verify the cells return to a single horizontal band — one row of five
   cells side by side (geometry round-trip) — with fewer rows and a
   different cell count than the two-row tiled grid shown just before.

### Scenario 2: Auto Layout gates control-panel visibility on viewer-box size

Steps:

1. Return to the Trellis Plot with DIS_POP as the single split column and
   Use Tiled View switched back on after Scenario 1's round-trip.
2. In the Properties panel, ensure Auto Layout is on (the default) and
   Show Control Panel is on.
3. Resize the Trellis Plot viewer box so that its height exceeds 300 px
   and its width exceeds 100 px; verify the control panel (the selector
   combo and property controls) is visible inside the viewer.
4. Shrink the viewer below the threshold (height <= 300 or width <= 100);
   verify the control panel is no longer shown inside the viewer. Confirm
   the Show Control Panel setting still reads its previously set value —
   the auto-hide must not change the stored setting, only what is
   displayed.
5. Restore the viewer box to a large size (height > 300, width > 100);
   verify the control panel is shown again without any explicit toggle.
6. Turn Auto Layout off in the Properties panel. With Show Control Panel
   explicitly on, shrink the viewer below the threshold; verify the
   control panel stays visible (explicit value wins at any size).
   Set Show Control Panel to off; verify the control panel is no longer
   shown, regardless of viewer size.
7. Re-enable Auto Layout and restore default settings for subsequent cleanup.

## Automation notes

- CHANNEL — `useTiledView`, `tilesPerRow`, `autoLayout`, `showControlPanel` and the split-column
  setup are prop writes, not Properties-panel typing (refdoc: Tiled View; a manual tester finds the
  first two in the Tiled View group of the gear Properties panel). The subject is the LAYOUT
  computed from those settings, so the assertion budget goes to geometry and visibility reads. UI
  actuation of the Tiled View group belongs to a separate ui-smoke pair.
- CHANNEL — Step 2 still clicks the gear (refdoc: Title bar): this section's only exercise of that
  hop, so a title-bar restructuring fails visibly instead of rotting.
- CHANNEL — the auto-hide threshold is driven with `page.setViewportSize()` (1920x1080 above,
  900x280 below), never by dragging the docking handle (refdoc: Hiding the control panel keeps the
  node).
- CHANNEL — `showControlPanel` is read off the property object as the other-channel witness that
  the size-driven auto-hide left the stored setting alone; the Context Panel would re-read the
  surface the auto-hide collapses.
- WITNESS — the control panel is reached through the `[name="viewer selector"]` parent hop and a
  "row" is inferred from cell top coordinates, both because the corresponding classes do not exist
  (refdoc: HTML Structure — No control-panel class, Cell rows). `maxPerRow` is the widest such
  group, bounded by both the category count and the viewport clamp.
- WITNESS — cell counts are `.d4-trellis-plot-cell` inside `[name="viewer-Trellis-plot"]`, every
  expectation derived from `col('DIS_POP').categories.length` through the viewport formulas.
- Y stays empty by design — a SINGLE split column is what engages `oneColumnOnly`.
- KNOWN GAP — the five-cell window shows 5 of the split column's 6 categories; the sixth needs a
  horizontal cell-grid scroll this section drives no channel for. A deliberate hole, not a `>= 5`
  hedge a regression to 5-of-any-width would also satisfy.

### Spec must keep

- The Scenario 1 ladder enters from the non-default side — Tiled View OFF before the Step 4 anchor,
  so Step 4 is a real off→on transition; else the on-entry write is a no-op and Step 4 witnesses
  nothing.
- The untiled entry is witnessed by a DOM read of the band's shape (one row, `min(5, cats)` cells)
  — else prop-echo.
- The control-panel witness stays TWO fields read in one pass off the `[name="viewer selector"]`
  parent — `exists` and computed `display` — asserted apart at every anchor (refdoc: Hiding the
  control panel keeps the node). Never fuse them into `!!cp && display !== 'none'` and never assert
  visibility by bare existence: a fused flag passes on the node-deletion regression the pair
  promises to catch.
- No `.d4-trellis-plot-control-panel` and no `.d4-trellis-plot-row` selectors (refdoc: HTML
  Structure) — such a locator matches nothing and passes a zero-count assert forever. Rows stay
  inferred from cell top coordinates, `showControlPanel` stays a property read.
- Step 7's untiled geometry stays a HORIZONTAL band — `rows` 1 plus `maxPerRow` equal to the
  rendered cell count (refdoc: The category viewport clamp). Never re-word it to a "vertical strip".
- Every expected cell count, row count and row width is a `Math.min(5, ...)` formula over the live
  category count (refdoc: The category viewport clamp) — never a literal 6, never the bare category
  count, or the expectation stops tracking the clamp this pair exists to catch.
- The split column stays WIDER than `tilesPerRow` 4, guarded by the live precondition assert
  (`expect(splitCats).toBeGreaterThan(DEFAULT_TILES_PER_ROW)`) rather than trusted by name — on a
  narrower column the two modes render the same geometry and Steps 4, 6 and 7 grade nothing.
- `count` / `maxPerRow` expectations are SLOT counts, a tiled grid padding its last row; the only
  category witness is Step 4's canvas-bearing-cell read at `min(cats, slots)`. "Fixing" a slot
  expectation down to the category count makes it false.
