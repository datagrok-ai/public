---
feature: matrixplot
realizes_atlas:
  - matrixplot.cp.axes-layout-and-style
  - matrixplot-auto-layout-overrides-axis-toggles
realizes:
  - viewers.matrix-plot
realized_as:
  - matrixplot-axes-layout-and-style-spec.ts
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19106
    status: fixed
  - id: GROK-18736
    status: fixed
expected_results:
  - anchor: "Show X Axes checkbox under Auto Layout"
    expectation: "With Auto Layout on (the default), unchecking the real Show X Axes
      checkbox collapses the shared X-axes tick strip canvas to 0x0 (the
      deterministic DOM signal), and re-checking it restores a strip canvas with
      a non-zero width (GROK-19106 guard, round-trip)"
  - anchor: "Explicit axes flags with Auto Layout off"
    expectation: "With Auto Layout off, unchecking Show X Axes collapses the top
      tick strip canvas to 0x0 and unchecking Show Y Axes collapses the left
      tick strip canvas to 0x0; re-checking each restores a non-zero strip
      canvas (the explicit flags act verbatim, round-trip)"
  - anchor: "Font property tracks into the label divs"
    expectation: "After changing the Font size from 10 to 24, the matrix
      column-label divs' computed font reads a 24px size (GROK-18736 guard);
      reverting the size to 10 returns the computed font to 10px"
---

# Matrix Plot — Axes Visibility, Auto Layout, Font

## Purpose

Verifies that the Show X Axes and Show Y Axes checkboxes really hide and
restore the shared axis tick strips — first under the default Auto Layout,
the exact configuration where GROK-19106 historically made the checkbox a
no-op, then with Auto Layout off, where the flags must act verbatim — and
that a Font property change visibly reaches the rendered column labels
(GROK-18736, the "style changes have no effect" class). Both checks read
deterministic DOM signals: the strip canvas size and the labels' computed
font; no pixel probing is needed.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Matrix plot to the current table view via the Toolbox viewer icon;
   the default 4x4 matrix on AGE, HEIGHT, WEIGHT, STARTED renders.

## Scenarios

### Scenario 1: Show X Axes checkbox with Auto Layout on (GROK-19106 guard)

Steps:
1. Open the viewer settings and verify **Auto Layout** (Style section) is on
   — the default.
2. Expand the **Axes** category in the settings panel — it is collapsed by
   default, so the axes rows are not visible until expanded.
3. Record the X-axes tick strip canvas size — it has a non-zero width while
   the axes are shown.
4. Click the real **Show X Axes** checkbox to uncheck it.
5. Verify the X-axes tick strip canvas collapsed to 0x0 — the deterministic
   DOM signal that the axes are hidden (the element stays in the DOM).
6. Click **Show X Axes** again to re-check it and verify the strip canvas is
   re-created with a non-zero width (round-trip).

Expected:
- With Auto Layout on, unchecking the Show X Axes checkbox collapses the top tick strip canvas to 0x0 and re-checking it restores a non-zero strip canvas

### Scenario 2: Explicit axes flags with Auto Layout off

Steps:
1. Turn **Auto Layout** (Style section) off.
2. Uncheck **Show X Axes** — the top tick strip canvas collapses to 0x0;
   re-check it — the strip canvas returns with a non-zero width.
3. Uncheck **Show Y Axes** — the left tick strip canvas collapses to 0x0;
   re-check it — the strip canvas returns with a non-zero width.
4. Turn **Auto Layout** back on (round-trip to the default).

Expected:
- With Auto Layout off, each Show X/Y Axes flag acts verbatim: unchecking collapses the matching tick strip canvas to 0x0 and re-checking restores it

### Scenario 3: Font change reaches the labels (GROK-18736 guard)

Steps:
1. Read the computed font of a matrix column-label div — the default size is
   10px.
2. In the settings **Font** row (Style section), change the size from 10 to
   24.
3. Verify the column-label divs' computed font now reads a 24px size — the
   style change visibly reached the rendered labels.
4. Revert the Font size to 10 and verify the computed font reads 10px again.

Expected:
- The label divs' computed font tracks the Font property: 24px after the change, 10px after the revert

## Automation notes

- GROK-19106 must be driven through the real checkbox, not a JS API property
  assignment — the property path never reproduced the bug. The Axes category
  is expanded by clicking the plus icon inside
  `[name="prop-category-axes"]`; the checkbox is
  `[name="prop-view-show-x-axes"]` (a real CDP click and a synthetic
  `.click()` both toggle it — use the real click to stay on the bug's path).
- The deterministic signal is
  `document.querySelector('[name="viewer-Matrix-plot"] .d4-layout-top canvas').width === 0`
  for the X strip and the `.d4-layout-left canvas` width for the Y strip. Do
  not run pixel probes on a collapsed strip — `getImageData` on the 0x0
  canvas throws IndexSizeError; check `canvas.width > 0` first.
- The Font row `[name="prop-view-font"]` is a composite editor: a size
  textbox plus a size list — a synthetic `.click()` on a size list item
  applies it, or type into the size textbox and dispatch input / Enter /
  change. The GROK-18736 signal is
  `getComputedStyle(labelDiv).font` on a column-label div (a leaf div whose
  trimmed text equals a column name), which tracks `props.font` — assert the
  computed style, not the property echo.
- The GROK-19106 recon was made on a large viewer; the Auto Layout
  size-heuristic can still hide axes on very small viewers regardless of the
  flags, so the spec keeps the viewer at its default docked size.

---
{
  "order": 13,
  "datasets": ["System:DemoFiles/demog.csv"]
}
