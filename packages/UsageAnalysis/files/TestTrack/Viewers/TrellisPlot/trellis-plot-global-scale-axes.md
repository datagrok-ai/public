---
feature: trellisplot
realizes_atlas:
  - trellisplot.cp.global-scale-inner-axes
  - trellisplot.int.global-scale-range-slider-sync
realizes:
  - viewers.trellis-plot
priority: p1
target_layer: playwright
boot_lane: local
coverage_type: regression
related_bugs:
  - id: GROK-14587
    status: fixed
realized_as:
  - trellis-plot-global-scale-axes-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After enabling Global Scale, at least two probed .d4-trellis-plot-cell
      canvases differ from their respective pre-toggle snapshots (settle-gated
      per-cell diff), and no console error is thrown.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      With Show Range Sliders on and the X or Y axis displayed, the right-click
      context menu on the trellis viewer contains the item 'Reset Inner Range
      Sliders'.
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      After hiding BOTH inner axes (Show X Axes and Show Y Axes off), the
      right-click context menu on the trellis viewer opens (it still carries
      'Properties...') and does NOT contain the item 'Reset Inner Range
      Sliders', confirming the item exists only while at least one inner-axis
      slider does.
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      After re-enabling both inner axes (Show X Axes and Show Y Axes on), the
      right-click context menu on the trellis viewer again contains 'Reset Inner
      Range Sliders'.
  - anchor: "Scenario 1 Step 11"
    expectation: >-
      The count of [type="range-slider"] elements under the trellis viewer root
      is strictly lower with an axis hidden than with both shown (a delta, not
      an absolute count).
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      After dragging the shared X-axis range slider to narrow the range, each
      probed .d4-trellis-plot-cell canvas differs from its pre-drag snapshot
      (settle-gated per-cell diff), and two distinct cells still differ from
      each other, confirming the shared re-bound applied to all cells; the drag
      throws no console error.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      After clicking 'Reset Inner Range Sliders' in the context menu, each
      probed .d4-trellis-plot-cell canvas returns EXACTLY to its pre-drag
      full-range baseline (settle-gated per-cell hash equal to the Step 1
      baseline, and therefore different from the narrowed state), and no console
      error is thrown.
  - anchor: "Scenario 3 Step 2"
    expectation: >-
      Allow Zoom reads off on the scatter plot inner viewer before anything is
      written to it (the untouched shipped default, read never repaired), and
      mouse-wheel scroll on the scatter plot cell produces no zoom change — the
      inner scatter plot axis bounds visible in the cell are identical before
      and after the wheel event (settle-gated per-cell canvas diff shows no
      change, and no console error is logged).
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      The switch to Bar chart lands — the viewer reports Bar chart as its inner
      type and the probed cell is repainted (its canvas differs from the scatter
      plot frame) — and mouse-wheel scroll on that bar chart cell then produces
      no zoom change: the cell canvas is identical before and after the wheel
      event (no zoom, no console error).
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      The switch to Box plot lands — the viewer reports Box plot as its inner
      type and the probed cell is repainted (its canvas differs from the bar
      chart frame) — and mouse-wheel scroll on that box plot cell then produces
      no zoom change: the cell canvas is identical before and after the wheel
      event (no zoom, no console error).
  - anchor: "Scenario 3 Step 5"
    expectation: >-
      After setting Allow Zoom to true on the scatter plot inner viewer via the
      inner-viewer property editor, mouse-wheel scroll inside a scatter plot
      cell produces a visible zoom change — the settle-gated per-cell canvas
      diff registers a non-trivial delta on the same cell before and after the
      wheel event, and no console error is thrown.
  - anchor: "Scenario 3 Step 6"
    expectation: >-
      After resetting Allow Zoom back to false (default) via the inner-viewer
      property editor, mouse-wheel scroll on the same scatter plot cell again
      produces no zoom change — the per-cell canvas diff shows no change and no
      console error.
---

# Trellis Plot — Global Scale, Inner Axes, and Range Slider Reset

## Setup

1. Close all open views.
2. Open System:DemoFiles/demog.csv.
3. Add a Trellis plot viewer (via the toolbar Add Viewer button or
   Toolbox > Viewers > Trellis plot).
4. In the Trellis plot control panel, set X to SEX and Y to RACE (2 × 4 = 8
   cells; all within the viewport clamp so cell count equals the category
   product).
5. In the control panel viewer-type selector, switch the inner viewer to
   Scatter plot.

## Scenarios

### Scenario 1: Global Scale toggle and axis-driven slider presence

Steps:
1. Note the current appearance of at least two trellis cells (for example,
   the cell in row 0 column 0 and the cell in row 1 column 0) to capture
   their current per-category axis bounds as a visual baseline.
2. In the Trellis plot Properties panel (gear icon or right-click > Properties),
   set Show X Axes and Show Y Axes to Always, so the inner axes are displayed
   regardless of the cell size the current split leaves (their default is Auto).
3. In Properties, set Show Range Sliders on. The scenario prescribes all three
   values explicitly rather than relying on defaults, so the slider-presence
   assertions below rest on a state the scenario itself established.
4. Observe that nothing has changed yet: the cells still show no inner axes and
   no range sliders. Neither axis setting nor Show Range Sliders has any effect
   while Global Scale is off — the next step turns it on.
5. Enable Global Scale in Properties (or via the control panel Global Scale
   checkbox if visible), then wait for the viewer to settle — all cells
   re-render with shared axis bounds.
6. Right-click the trellis viewer background (not inside a cell) to open its
   context menu.
7. Verify that the context menu contains the item 'Reset Inner Range Sliders'.
8. In Properties, set BOTH Show X Axes and Show Y Axes to off. Wait for the
   viewer to settle, then right-click the trellis viewer background again.
9. Confirm the context menu did open — it still lists 'Properties...' — and then
   verify that 'Reset Inner Range Sliders' is NOT present in it.
10. Set both Show X Axes and Show Y Axes back to on, wait for the viewer to
    settle, right-click the trellis viewer background once more, and verify that
    'Reset Inner Range Sliders' is present again.
11. With both axes shown, count the range slider controls visible inside the
    trellis viewer and record this as the axes-shown count. Then set Show X Axes
    to off, wait for the viewer to settle, and count again — the axes-hidden
    count. Verify the axes-shown count is greater than the axes-hidden count.
    Set Show X Axes back to on before continuing to Scenario 2.

### Scenario 2: Shared range slider re-bounds all cells; Reset restores full range

Steps:
1. Starting from the state at the end of Scenario 1 (Global Scale on, Show
   Range Sliders on, Show X Axes on, Scatter plot inner viewer, X = SEX,
   Y = RACE). Note the current appearance of the same two trellis cells used
   in Scenario 1 as their full-range baseline.
2. Locate the shared X-axis range slider rendered between the inner-viewer
   cells (it appears as a horizontal slider at the viewer level, separate from
   the category scroll controls on the viewer edges).
3. Verify that the two observed cells currently look different from each other
   (different category subsets should produce different scatter distributions).
4. Drag the shared X-axis range slider to narrow it to approximately the
   central half of the full data range for the X column.
   Wait for the viewer to settle.
5. Right-click the trellis viewer background, select 'Reset Inner Range Sliders'.
   Wait for the viewer to settle.

### Scenario 3: Inner viewers do not zoom on mouse-wheel (GROK-14587)

Regression guard for GROK-14587: inner viewers (scatter plot, bar chart, box
plot) must not zoom on mouse-wheel by default. The mechanism differs by inner
viewer type, and the scenario says so rather than pretending it is uniform:

- **Scatter plot has a gate.** An Allow Zoom property exists on the scatter
  look; it is on for a standalone scatter plot and the trellis preset turns it
  off, which is exactly the surface GROK-14587 guards. Its opt-in round-trip
  (Steps 5-6) is a scatter-plot-only contract.
- **Bar chart and box plot have no gate at all.** Neither exposes an Allow Zoom
  property — their inner-viewer property tab has no such row — and neither
  zooms on the wheel. Steps 3-4 therefore verify the absence of zoom itself, not
  the state of a setting.

Coverage stays on all three inner viewer types: the ticket is about all three
not zooming, and the wheel probe is the same for each.

Steps:
1. Bring the viewer to this scenario's own start state: in the Trellis plot
   Properties panel turn Global Scale off and Show Range Sliders off, keep
   X = SEX and Y = RACE, and set the inner viewer type to Scatter plot in the
   control panel viewer-type selector. Wait for all 8 cells to render.
   (The inner range sliders appear only while the axis area is hovered, and
   every wheel gesture below has to hover a cell — leaving them on would move
   pixels that the per-cell comparison is watching.)
2. Open the inner-viewer property editor (gear icon on the trellis panel title
   bar, then the Scatter plot tab, expanding the collapsed Misc category) and
   READ the Allow Zoom value without changing it — it must already be off. Do
   not click the checkbox: setting it would repair the very default this step
   checks. Note the current appearance of the top-left cell
   (SEX = F, RACE = Asian), hover the mouse over the centre of that cell and
   scroll the mouse wheel downward five steps, then wait for the viewer to
   settle.
3. In the control panel switch the inner viewer type to Bar chart and wait for
   the cells to settle. Confirm the switch really happened before probing: the
   viewer now names Bar chart as its inner type and the top-left cell is
   repainted — it no longer shows the scatter plot it showed in Step 2. There is
   no Allow Zoom row to look at for this inner type. Note the top-left cell's
   appearance, then scroll the mouse wheel downward five steps over its centre
   and wait for the viewer to settle.
4. In the control panel switch the inner viewer type to Box plot and wait for
   the cells to settle. Again confirm the switch landed: the viewer names Box
   plot as its inner type and the top-left cell is repainted from the bar chart
   it showed in Step 3. As with Bar chart, there is no Allow Zoom row for this
   inner type. Note the top-left cell's appearance, then scroll the mouse wheel
   downward five steps over its centre and wait for the viewer to settle.
5. Switch the inner viewer type back to Scatter plot. In the inner-viewer
   property editor set Allow Zoom to true and wait for the cells to settle.
   Note the top-left cell's appearance, then scroll the mouse wheel downward
   five steps over its centre and wait for the zoom to settle.
6. In the inner-viewer property editor set Allow Zoom back to false and wait
   for the cells to settle. Note the top-left cell's appearance again, then
   scroll the mouse wheel downward five steps over its centre and wait for the
   viewer to settle.

## Automation notes

- CHANNEL — Allow Zoom is driven through the real property-grid editor
  (`tr[name="prop-allow-zoom"]`), not props; the spec expands the owning category header first,
  guarded against a re-entry collapsing it (refdoc: pitfall 16).
- CHANNEL — the inner type switch is a setting, not the tested action, so it goes through
  `props.viewerType`; only the wheel is driven as real input.
- CHANNEL — the shared-slider drag goes through the trusted-input recipe (refdoc: Dragging an
  inner-axis slider), never a synthetic sequence.
- CHANNEL — every context-menu label read is scoped to `.d4-menu-popup .d4-menu-item-label`
  (refdoc: pitfall 22); an unscoped sweep grades a different menu.
- CHANNEL — read and write are separate helpers: `openInnerViewerTab` is navigation glue and writes
  nothing, `allowZoomState` only reads, `setAllowZoom` is reserved for the explicit Steps 5-6
  transitions. Its gear-open third is the registered helper `openViewerGear` (`helpers/viewers.ts`);
  only the tab-select and Misc-expand halves are spec-local.
- CHANNEL — no-error floor: `pageerror` + `console(error)` subscribed before login, benign ambient
  classes filtered, each step asserting a DELTA across its own actuation. These two are the only
  channels available (refdoc: pitfall 26).
- target_layer playwright: context-menu presence, slider element counts, a slider drag and
  settle-gated per-cell canvas diffs are all unreachable from the apitest layer.
- WITNESS — inner-axis slider assertions are a COUNT of `[type="range-slider"]` scoped to the
  trellis viewer root (refdoc: Discriminating inner-axis sliders), graded as a shown-vs-hidden
  DELTA: that selector also matches the category scroll sliders, so an absolute count or a
  presence check grades the wrong elements.
- WITNESS — the `Reset Inner Range Sliders` absence check hides BOTH axes; the body, the anchors
  and the spec all encode the both-axes form, because the item's gate is an OR over the two sliders
  (refdoc: Discriminating inner-axis sliders) and a one-axis hide would grade nothing.
- WITNESS — the identity half of each inner-type switch is a `props.viewerType` read-back, weak
  alone as a same-channel echo, paired with a per-cell canvas delta: these three inner types put no
  named node in the cell to look for (refdoc: A cell hosts the inner viewer's bare CANVAS). The
  delta is the render-signal index's prescribed signal for the props path, which never fires
  `d4-trellis-plot-viewer-type-changed`.
- WITNESS — the Bar chart and Box plot arms read NO property: `allowZoom` does not exist on those
  looks (refdoc: Allow Zoom), so their only signal is the canvas staying byte-identical across the
  wheel.
- Column prescriptions SEX and RACE (refdoc: pitfall 19) keep the 8-cell grid under the per-axis
  clamp (refdoc: The category viewport clamp), so the DOM cell count always equals the category
  product.
- Scenario 3 Step 1 resets `globalScale` and `showRangeSliders` to false: the Scenario 1-2 ladder
  leaves them on, the sliders are hover-revealed, and every wheel gesture must hover a cell first.
- Scenario 3 carries the GROK-14587 guard; the bug id lives in `related_bugs`, not
  `realizes_atlas`.
- `int.global-scale-range-slider-sync` (one shared slider re-bounds all cells, Reset restores
  every cell) is closed by Scenario 2.

### Spec must keep

- Null-guard BOTH points of every canvas delta — else the Step 5 "changed" assert passes on a
  vanished canvas and the Steps 2/3/4/6 "unchanged" asserts compare null to null.
- The wheel stays on trusted `page.mouse.move` + `page.mouse.wheel` — synthetic wheel events never
  reach the d4 handler, making every "unchanged" assertion vacuous.
- The Allow Zoom round-trip true → false is mandatory: Step 5 alone is only a positive control,
  Step 6 is what proves the gate rather than the wheel plumbing changed.
- All three inner types (Scatter plot, Bar chart, Box plot) stay in the no-zoom sweep — the wheel
  zoomed ACROSS inner types, so Scatter alone cannot cover it.
- Every inner-type switch in Scenario 3 Steps 3-4 keeps a hard landing witness: the type read back
  PLUS a per-cell canvas delta — else a no-op switch leaves both arms re-probing the SCATTER cell,
  both `after === before` asserts pass and the Bar/Box coverage vanishes with no red step.
- The cell count is NOT that witness — 8 is the category product for every inner type, so
  `toHaveCount(8)` is a settle wait. A DOM node lookup is no substitute either: no named viewer
  node exists for these types.
- The repaint delta may never be dropped: it is the half independent of the channel the switch was
  written through.
- Scenario 2 Step 5 grades Reset as an EXACT per-cell return to the pre-drag baseline (refdoc:
  Discriminating inner-axis sliders), never as a move away from the narrowed state — the weak form
  is satisfied by any repaint the context-menu dispatch happens to cause.
- Read context-menu labels ONLY through `.d4-menu-popup .d4-menu-item-label` (refdoc: pitfall 22)
  — unscoped, the Step 9 witness `Properties...` rests on that string happening to be unique
  document-wide.
- The Step 9 absence check stays gated on a positive witness that the menu opened (`Properties...`,
  present regardless of axis state) — an unopened menu yields an empty label list, and an empty
  list satisfies `not.toContain`.
- Step 2 reads the Allow Zoom default UNTOUCHED and asserts it before any write; a setter first,
  even a "no-op if already false" one, repairs a regression instead of catching it. Navigation to
  the row is allowed, writing is not.
- Step 2 grades the untouched default strictly `toBe(false)`, never `toBeFalsy` / `!x` /
  `not.toBeChecked` — `allowZoomState` returns null on a missing row, and null passes those forms.
- Do not add an Allow Zoom read to the Bar chart or Box plot arms: the property does not exist
  there, so the assertion is either vacuous or aimed at a missing row.
- The console floor stays a per-step delta over the `pageerror` / `console(error)` collectors — a
  global `toEqual([])` is false-red on noise accumulated by earlier steps.
- No assertion may read `grok.shell.warnings`: the symbol does not exist, so
  `(grok.shell.warnings ?? []).length` is 0 on both sides and can never fail. Do not reintroduce it
  as a "secondary proxy".
