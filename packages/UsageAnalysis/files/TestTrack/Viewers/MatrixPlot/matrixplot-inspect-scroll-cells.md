---
feature: matrixplot
realizes_atlas:
  - matrixplot.cp.inspect-scroll-cells
  - matrixplot-large-matrix-scroll-viewport
  - matrixplot-cell-tooltip-then-open-fullscreen
realizes:
  - viewers.matrix-plot
realized_as:
  - matrixplot-inspect-scroll-cells-spec.ts
priority: p1
target_layer: playwright
coverage_type: smoke
related_bugs: []
expected_results:
  - anchor: "Slider drag moves the viewport"
    expectation: "Dragging the horizontal viewport slider's max handle toward the
      slider midpoint narrows the visible X columns (fewer inner cell canvases
      and fewer top label divs than the full 4; the exact settle count is
      viewport-dependent); dragging the handle back restores 16 cells and the
      full four X labels (round-trip)"
  - anchor: "Viewport cell cap rejection"
    expectation: "On the 16x16 column set (the four base demog columns plus 12
      fixture-derived numeric columns), opening both viewport sliders fully is
      rejected at the 250-cell cap: every sampled visible cell count stays <=
      250 and never reaches 256; the exact settle value is viewport-dependent
      (224 and 240 observed live) — the over-cap drag increments leave the
      viewport unchanged"
  - anchor: "Cell hover tooltip and expand icon"
    expectation: "Hovering an off-diagonal cell populates the tooltip element with
      the cell's identity text naming its X and Y columns (read only while the
      tooltip rectangle has a non-zero width), and the hovered cell's expand
      icon becomes visible (its visibility:hidden style is removed)"
  - anchor: "Expand opens a matching standalone viewer"
    expectation: "Clicking the expand icon on an off-diagonal cell with Cell Plot
      Type Density plot adds a standalone Density plot viewer to the view
      (viewer-set membership), and clicking the expand icon on a diagonal cell
      adds a standalone Histogram; both added viewers are closed afterwards"
  - anchor: "Per-cell wheel zoom"
    expectation: "A wheel zoom over one cell changes that cell's settle-gated pixel
      measurement while a neighbouring cell's measurement stays unchanged, and
      the zoom raises no new console or page error"
---

# Matrix Plot — Viewport Scrolling and Cell Inspection

## Purpose

Verifies how a user explores the matrix: moving the visible window with the
viewport sliders (and, on a fixture-extended 16x16 matrix, the 250-cell cap
rejecting an over-cap drag), hovering
a cell to read which column pair it plots and to reveal its expand icon,
opening a cell as a standalone full-size viewer of the matching type, and
zooming inside a single cell. The tooltip text, the expand icon's
visibility, and the added viewer's type are readable product-state signals;
the zoom is a per-cell canvas change measured against an untouched
neighbouring cell.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Derive 12 numeric fixture columns on demog via the JS API (float columns
   computed as AGE + i, i = 1..12) so the eligible numeric/datetime set
   reaches 16 columns. These fixture columns are removed in a `finally`
   teardown at the end of the run.
4. Add a Matrix plot to the current table view via the Toolbox viewer icon
   and set its X and Y column sets to the four base columns AGE, HEIGHT,
   WEIGHT, STARTED (via props); the 4x4 matrix renders with
   **Cell Plot Type** Density plot.

## Scenarios

### Scenario 1: Move the viewport with the sliders

Steps:
1. Verify the horizontal and vertical viewport sliders are present on the
   4x4 matrix.
2. Drag the horizontal slider's max handle toward the slider midpoint.
3. Verify the viewport moved: the visible X columns narrowed — fewer than 4
   inner cell columns remain and the top label-div set shrinks below the four
   full labels (the exact settle count is viewport-dependent).
4. Drag the max handle back to the slider end and verify the full 4x4
   viewport (16 cells, four X labels) is restored (round-trip).
5. Build the 16x16 matrix: set both the X and Y column sets to all 16
   numeric columns (the four base plus the 12 fixture columns) via props.
   The re-tiled viewport initially shows only a 5x5 window (25 cells).
   Drag the horizontal slider's max handle fully open, then the vertical
   slider's max handle fully open, trying to expose all 256 cells. Verify
   the over-cap drag is REJECTED at the 250-cell cap: every sampled visible
   cell count stays <= 250 and never reaches 256; the exact settle value is
   viewport-dependent (224 and 240 observed live) — the viewport stays put
   on the final over-cap increments.

Expected:
- The slider drag narrows the visible X columns to fewer than 4 (at least 1; the exact count is viewport-dependent) and back, and on the 16x16 set the over-cap drag is rejected: every sampled visible cell count stays <= 250 and never reaches 256

### Scenario 2: Hover a cell — tooltip and expand icon

Steps:
1. Hover an off-diagonal cell (for example X=HEIGHT, Y=AGE).
2. Verify the tooltip appears and its text names the hovered cell's X and Y
   columns (X: HEIGHT, Y: AGE). Read the tooltip text only while the
   tooltip rectangle has a non-zero width.
3. Verify the hovered cell's expand icon became visible — before the hover
   it is present in the DOM but hidden.

Expected:
- The tooltip identifies the hovered cell's X and Y columns and the cell's expand icon becomes visible on hover

### Scenario 3: Open a cell as a standalone viewer

Steps:
1. Record the view's current viewer set.
2. Hover an off-diagonal cell and click its expand icon.
3. Verify a standalone Density plot viewer was ADDED to the view — the added
   viewer's type matches the configured **Cell Plot Type** (viewer-set
   membership, not a pixel check).
4. Close the added viewer.
5. Hover a diagonal cell (for example X=AGE, Y=AGE) and click its expand
   icon.
6. Verify a standalone Histogram viewer was added — diagonal cells always
   plot a histogram.
7. Close the added viewer and verify the viewer set returned to the recorded
   state.

Expected:
- The off-diagonal expand adds a Density plot (matching Cell Plot Type), the diagonal expand adds a Histogram, and closing both restores the original viewer set

### Scenario 4: Zoom inside one cell

Steps:
1. Note the current console-error and page-error counts. Measure settle-gated
   pixel readings of two off-diagonal cells: the zoom target and an untouched
   neighbour.
2. Send a wheel zoom-in over the target cell's center; wait for the redraw.
3. Verify the target cell's pixel reading changed while the neighbour's
   reading stayed unchanged — the zoom is per-cell.
4. Send a wheel zoom-out over the same cell to move back toward the original
   scale.
5. Verify no new console or page error appeared across the zoom gestures.

Expected:
- The wheel zoom changes only the target cell's settle-gated pixel measurement (the neighbour is the control group) and raises no new console or page error

## Automation notes

- Slider handles: `svg[name="x-slider"]` / `svg[name="y-slider"]` with
  children `[name="pan-handle"]`, `[name="min-handle"]`, `[name="max-handle"]`.
  Synthetic dispatched mousedown -> document-mousemove steps ->
  document-mouseup DO drive the drag (same widget family as the PC Plot axis
  sliders). The sliders are present already at 4x4.
- The 250-cell cap portion runs on a fixture-extended demog: 12 derived
  float columns (AGE + i, i = 1..12) are added via the JS API in Setup so
  the eligible numeric set reaches 16. Recon on the 16x16 set: the initial
  re-tiled viewport renders a small window (fewer than 80 cells; 25 observed);
  opening the x slider fully grows the count (viewport-dependent); opening
  the y slider fully settles below the cap, NEVER 256 — the cap rejects the
  final increments and the drag stops short. Assert the invariant: every
  sampled visible cell count stays <= 250 and never reaches 256; the exact
  settle value is viewport-dependent (224 and 240 observed live). The fixture
  columns are removed from the dataframe in a `finally` teardown so the
  extended demog never leaks past the run.
- Cell indexing is row-major over the VISIBLE sets:
  `cells[i]` maps to `(y = floor(i / nX), x = i mod nX)`; with identical X
  and Y sets the diagonal cells are those with `floor(i / nX) === i mod nX`.
- The hover is synthetic `mouseenter` + `mousemove` at the cell-canvas
  center. The tooltip is the document-level `.d4-tooltip`; it RETAINS its
  last text when hidden, so always gate the read on
  `getBoundingClientRect().width > 0`. The expand icon
  (`[name="icon-expand-arrows"]` inside the cell's parent) is ALWAYS in the
  DOM with `visibility:hidden` until hover — assert `style.visibility`, not
  presence (`offsetParent` stays truthy while hidden).
- The expand click is a synthetic `.click()` on the revealed icon; the
  membership signal is `grok.shell.tv.viewers.map(v => v.type)` gaining
  'Density plot' (or 'Scatter plot' when that type is configured) /
  'Histogram'. The added viewer is closed via its `close()`.
- The wheel zoom is a synthetic `WheelEvent` (`deltaY: -300` in, `+300`
  out) at the cell-canvas center; the zoom-out is approximate, so only the
  zoom-in delta and the neighbour's invariance are asserted, plus the
  no-error floor. Cell canvas backing-store pixels differ from CSS px —
  scale coordinates by `rect.width / canvas.width`.
- In-cell POINT interactions (hover-highlight, click-to-current-row,
  drag-select) are deliberately NOT covered anywhere for this viewer — an
  operator decision; live recon also found them non-functional via trusted
  input on the reference build. Do not add selection asserts against inner
  cells.

---
{
  "order": 15,
  "datasets": ["System:DemoFiles/demog.csv"]
}
