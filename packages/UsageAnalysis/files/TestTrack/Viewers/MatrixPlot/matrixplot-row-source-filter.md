---
feature: matrixplot
realizes_atlas:
  - matrixplot.cp.row-source-filter
  - matrixplot-row-source-filter-narrow-cells
realizes:
  - viewers.matrix-plot
realized_as:
  - matrixplot-row-source-filter-spec.ts
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs: []
expected_results:
  - anchor: "Row Source Selected repaints over the selection"
    expectation: "With exactly 50 rows selected (precondition df.selection.trueCount
      == 50), switching Row Source to Selected changes the matrix's settle-gated
      pixel measurement versus the Filtered baseline; switching Row Source back
      to Filtered and clearing the selection returns the measurement EXACTLY to
      the baseline (round-trip)"
  - anchor: "Filter Panel narrows the shared filter"
    expectation: "Applying a SEX = M categorical filter on the Filter Panel drops
      df.filter.trueCount below the full row count; removing the filter restores
      df.filter.trueCount to the full count (a product-state round-trip)"
  - anchor: "Viewer Filter property repaints without touching the shared filter"
    expectation: "Setting the viewer's Filter property to `${AGE} > 30` changes the
      matrix's settle-gated pixel measurement while df.filter.trueCount is NOT
      asserted for this step (the viewer-level formula does not drive the shared
      filter); clearing the Filter property restores the pixel measurement
      EXACTLY to the baseline"
---

# Matrix Plot — Row Source and Filtering

## Purpose

Verifies that narrowing the plotted rows reaches every inner cell of the
matrix through its three routes: the Row Source property switched to the
current selection, a Filter Panel filter on the shared table, and the
viewer's own Filter formula. The inner cells draw to canvas, so the redraws
are measured as settle-gated pixel changes with an exact restore after each
round-trip; the Filter Panel route also has a product-state signal — the
shared table's filtered row count — which is asserted only for that route,
because the viewer-level formula deliberately does not touch the shared
filter.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Matrix plot to the current table view via the Toolbox viewer icon;
   the default 4x4 matrix renders.

## Scenarios

### Scenario 1: Row Source switched to Selected

Steps:
1. Measure a settle-gated pixel reading of the matrix in the default state
   (**Row Source** Filtered, empty selection) — the baseline.
2. Select the first 50 rows of the table and verify the precondition: the
   selected row count equals exactly 50.
3. Set **Row Source** (Data section) to Selected; wait for the redraw.
4. Verify the pixel reading changed versus the baseline — the matrix redrew
   over only the selected rows.
5. Set **Row Source** back to Filtered and clear the selection; wait for the
   redraw.
6. Verify the pixel reading returned EXACTLY to the baseline (round-trip —
   rendering is deterministic for identical state at a fixed viewer size).

Expected:
- With 50 rows selected, Row Source Selected changes the settle-gated pixel measurement; reverting to Filtered and clearing the selection restores it exactly

### Scenario 2: Filter Panel filter narrows the shared table

Steps:
1. Record the full filtered row count of the table.
2. Open the Filter Panel and apply a categorical filter on **SEX**, keeping
   only **M**.
3. Verify the filtered row count dropped below the full count — the panel
   filter drives the shared table filter the matrix listens to.
4. Remove the SEX filter from the Filter Panel.
5. Verify the filtered row count restored to the full count (round-trip).

Expected:
- The SEX = M panel filter drops the filtered row count below the full count; removing it restores the full count

### Scenario 3: Viewer Filter property

Steps:
1. Measure a settle-gated pixel reading of the matrix with the **Filter**
   property (Data section) empty — the baseline.
2. Set the **Filter** property to `${AGE} > 30`; wait for the redraw.
3. Verify the pixel reading changed versus the baseline. Do NOT read the
   shared table's filtered row count here — the viewer-level formula narrows
   only this viewer's rows and leaves the shared filter untouched.
4. Clear the **Filter** property; wait for the redraw.
5. Verify the pixel reading returned EXACTLY to the baseline (round-trip).

Expected:
- The `${AGE} > 30` Filter property changes the settle-gated pixel measurement and clearing it restores the measurement exactly; the shared filtered row count is not asserted on this route

## Automation notes

- The viewer handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Matrix plot')`; row counts
  are read dynamically from `df.filter.trueCount` and
  `df.selection.trueCount`, never hard-coded.
- The pixel measurement is a settle-gated hash / non-white pixel count over
  the inner cell canvases
  (`[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer`),
  re-read every ~300 ms until two consecutive readings agree (settle is
  ~1.5–2 s on demog). Exact-restore asserts are valid because rendering is
  deterministic for identical state at a fixed viewer size — keep the viewer
  size unchanged across a round-trip.
- Scenario 1 drives the selection via the JS API
  (`df.selection.init(i => i < 50)`) and Row Source via
  `mp.props.rowSource = 'Selected'` / `'Filtered'`; the panel row is
  `[name="prop-view-row-source"]` when a UI drive is preferred.
- Scenario 2 applies the panel filter via
  `grok.shell.tv.getFiltersGroup().updateOrAdd({type: 'categorical',
  column: 'SEX', selected: ['M']})` and removes it through the same group;
  the signal is `df.filter.trueCount`, asserted only on this Filter-Panel
  route.
- Scenario 3 sets `mp.props.filter = '${AGE} > 30'` and clears it with
  `''`. Live recon confirmed the viewer filter prop does NOT change
  `df.filter.trueCount` — asserting it here would be wrong by design, so the
  route is measured by pixels only.
- Spec floor calibration (whole-matrix ink deltas; re-calibrate if the viewer
  size changes): Row Source Filtered -> Selected (50 of 5850 rows) — delta
  ~17000, spec floor 2000; viewer Filter `${AGE} > 30` — delta ~600, spec
  floor 100.

---
{
  "order": 14,
  "datasets": ["System:DemoFiles/demog.csv"]
}
