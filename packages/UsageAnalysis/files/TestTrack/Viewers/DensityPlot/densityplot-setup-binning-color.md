---
feature: densityplot
realizes_atlas:
  - densityplot.cp.setup-binning-color
realizes:
  - viewers.density-plot
priority: p0
target_layer: playwright
coverage_type: smoke
realized_as:
  - densityplot-setup-binning-color-spec.ts
related_bugs:
  - id: GROK-16612
    status: fixed
  - id: GROK-17118
    status: fixed
expected_results:
  - anchor: "Add the viewer and pick the axis columns through the on-viewer selectors"
    expectation: "After picking a column through the on-viewer X or Y selector
      popup, the viewer's xColumnName / yColumnName props read the picked column
      (a product read of the UI actuation, not a prop echo) and the selector's
      DOM label shows the same name; no new browser console errors appear across
      several axis switches (GROK-16612)"
  - anchor: "Bin count and bin shape"
    expectation: "Changing Bins from 5 to 200 produces a strong settle-gated canvas
      repaint — the non-white pixel fraction shifts by a clear live-calibrated
      margin (reference ~0.43 at bins=5 vs ~0.08 at bins=200); the
      rectangle/hexagon bin-shape switch repaints or, if the live delta is too
      modest, holds an error-free floor"
  - anchor: "Bin-count lower boundary — bins = 1"
    expectation: "Setting Bins to 1 (the slider's lower clamp) collapses the plot to
      a single filled cell, so the non-white pixel fraction rises sharply above
      the 50-bin grid by a live-calibrated margin (reference px1 ~300972 vs px50
      ~158509, delta ~142463) with no page or console errors; reverting to 50
      bins drops the fill back (round-trip)"
  - anchor: "Invert Color Scheme and Color Transform Type"
    expectation: "Turning Invert Color Scheme on produces a strong settle-gated
      canvas repaint (reference non-white fraction 0.27 -> 0.81 with rectangle
      bins); switching Color Transform Type from linear to logarithmic repaints
      by a live-calibrated luminance delta or holds an error-free floor"
  - anchor: "GROK-17118 guard — the X selector offers numeric columns only"
    expectation: "With the rectangle bin shape active, searching a non-numeric
      column (SEX) in the X selector popup matches nothing and Enter commits
      nothing — the X column stays unchanged and no page or console errors
      appear (the numeric-only selector IS the fixed behaviour; the JS-API
      invalid-column path is deliberately not driven)"
  - anchor: "Degenerate zero-width bin range"
    expectation: "Assigning a numeric constant column (min == max, so a zero-width
      bin range) to the X axis through the on-viewer selector and rendering both
      rectangle and hexagon bin shapes completes with no page or console errors
      and no freeze — the viewer stays attached and each shape produces a
      settled non-zero repaint (a degenerate strip; the divide-by-zero bin-size
      path is guarded); the fixture column is removed afterwards even on
      failure"
  - anchor: "Layout and project persistence at peak configuration — layout round-trip"
    expectation: "After saving the layout of the fully configured viewer (X=AGE,
      Y=HEIGHT, bins=200, rectangle shape, inverted color scheme), adding a
      Scatter plot, and re-applying the saved layout, the view's viewer set
      equals the SAVED set (a Density Plot is present AND the later-added
      Scatter plot is absent) AND the restored viewer carries xColumnName AGE,
      yColumnName HEIGHT, bins 200, binShape rectangle, invertColorScheme true;
      the probe layout is deleted afterwards even on failure"
  - anchor: "Layout and project persistence at peak configuration — project save /
      Close All / reopen"
    expectation: "After saving the view as a project through the real ribbon Save
      button, Close All, and reopening the saved project, a Density Plot viewer
      is present in the reopened view AND its xColumnName, yColumnName, bins,
      binShape and invertColorScheme equal the persisted values (a cross-session
      round-trip); the probe project is deleted afterwards even on failure"
---

# Density Plot — Setup, Axis Columns, Binning, Color Mapping, Persistence

## Purpose

Verifies the Density Plot's day-to-day configuration surface on the demog
dataset: adding the viewer, picking the X and Y columns through the on-viewer
selectors (the real UI path), changing the bin count and bin shape, inverting
the color scheme and switching the color transform, the numeric-only column
guard, and the survival of a fully configured plot across a saved layout and a
saved project. The plot draws everything to a single canvas, so rendered
outcomes are checked as settle-gated pixel deltas with live-calibrated
thresholds, or as an honest error-free floor where the delta is too modest to
trust.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.

## Scenarios

### Scenario 1: Add the viewer and pick the axis columns through the on-viewer selectors

Steps:
1. In the Toolbox, click the **Density Plot** icon — the viewer attaches to the
   demog table view (on demog it auto-picks X=AGE, Y=HEIGHT).
2. Note the current page-error and console-error counts.
3. Click the on-viewer **X column selector** (bottom center of the plot) — a
   column popup opens; type WEIGHT, press ArrowDown, then Enter.
4. Verify the viewer's X column reads WEIGHT (product read of `xColumnName`)
   and the X selector label shows WEIGHT.
5. Click the on-viewer **Y column selector** (vertical, left edge) and pick
   WEIGHT the same way; then switch Y back to HEIGHT.
6. Switch X back to AGE through the X selector.
7. Verify the final state: X column AGE, Y column HEIGHT — both the props and
   the selector labels agree.
8. Verify no new page or console errors appeared across all the axis switches
   (GROK-16612 guard: switching axis columns via the on-viewer selectors used
   to error to the console on every switch).

Expected:
- Each pick through an on-viewer selector updates the matching column prop and the selector's DOM label
- No new browser console errors across several axis switches (GROK-16612)

### Scenario 2: Bin count and bin shape

Steps:
1. Open the viewer's settings (gear icon) — the property panel opens.
2. Bin-count lower boundary: set **Bins** to 50 (the default) and measure the
   canvas ink, then set **Bins** to 1 — the slider's lower clamp collapses the
   plot to a single filled cell; measure again.
3. Verify the bins=1 ink rises sharply above the 50-bin grid by a
   live-calibrated margin, with no page or console errors.
4. Revert **Bins** to 50 and verify the fill drops back (round-trip).
5. Set **Bins** to 5 — the plot repaints with a few large bins; measure the
   canvas ink.
6. Set **Bins** to 200 — the plot repaints with a fine grid; measure again.
7. Verify a strong render change between the two states (the bins 5 vs 200
   repaint is the strongest pixel signal on this viewer).
8. Set **Bin Shape** to rectangle (the default is hexagon) — the bins change
   form.
9. Verify a render change from the shape switch, or that the switch completes
   with no page or console errors if the live-calibrated delta is too modest.
10. Leave the configuration as is: bins 200, rectangle shape (nothing is
   reverted — the persistence tail persists this peak configuration).

Expected:
- Bins 5 -> 200 produces a strong settle-gated canvas repaint (live-calibrated threshold)
- Bins 1 (the lower clamp) inks sharply above the 50-bin grid with no errors, and reverting to 50 drops it back
- The rectangle/hexagon switch repaints or holds an error-free floor

### Scenario 3: Invert Color Scheme and Color Transform Type

Steps:
1. In the **Style** section of the property panel, turn **Invert Color
   Scheme** on — the color mapping of the binned area inverts; measure the
   canvas ink before and after.
2. Verify a strong render change (with rectangle bins the inversion flips the
   background of the binned area).
3. Set **Color Transform Type** to logarithmic — the count-to-color mapping
   changes.
4. Switch **Color Transform Type** back to linear, then to logarithmic again.
5. Verify a render change from the transform switches, or that the sequence
   completes with no page or console errors if the live-calibrated delta is
   too modest.
6. Leave Invert Color Scheme on and Color Transform Type logarithmic (nothing
   is reverted before the persistence tail).

Expected:
- Invert Color Scheme produces a strong settle-gated canvas repaint
- The linear/logarithmic transform switches repaint or hold an error-free floor

### Scenario 4: GROK-17118 guard — the X selector offers numeric columns only

The rectangle bin shape used to error when a column of an invalid
(non-numeric) type was assigned. The fixed behaviour is that the on-viewer
column selectors offer numeric columns only, so the invalid assignment cannot
be made through the UI. This guard drives the real UI controls only — the
JS-API invalid-column path is deliberately out of scope.

Steps:
1. Confirm the **Bin Shape** is rectangle (set in Scenario 2).
2. Note the current page-error and console-error counts.
3. Click the on-viewer **X column selector** — the column popup opens.
4. Type SEX (a non-numeric demog column) into the popup search.
5. Verify no matching column is offered; press Enter and verify nothing is
   committed — the popup stays open and the X column is unchanged.
6. Close the popup by clicking outside it.
7. Verify the X column still reads AGE and no new page or console errors
   appeared.

Expected:
- The X selector popup offers no non-numeric column: the SEX search matches nothing, Enter commits nothing, and the X column stays AGE with no errors

### Scenario 5: Degenerate zero-width bin range

A numeric constant column has an identical minimum and maximum, so its bin
range has zero width. The bin-size divisor is guarded against divide-by-zero;
this edge confirms the guard holds through the real UI for both bin shapes.
The constant column is a fixture created via the JS API and removed in a
teardown that runs even when a verification fails.

Steps:
1. Create a numeric constant fixture column (every value 42) on the demog
   table through the JS API.
2. Note the current page-error and console-error counts.
3. Assign the constant column to the X axis through the on-viewer **X column
   selector** — a numeric column is offered, so the pick commits.
4. Set **Bin Shape** to rectangle and verify a settled repaint completes: the
   viewer stays attached and paints a non-zero degenerate strip with no page
   or console errors and no freeze.
5. Set **Bin Shape** to hexagon and verify the same error-free floor.
6. Restore the X axis to AGE and the bin shape to rectangle, then remove the
   fixture column (teardown runs even on failure).

Expected:
- Assigning a zero-width (constant) column to X and rendering both bin shapes completes with no errors and no freeze — the viewer stays attached and each shape produces a settled non-zero repaint; the fixture column is removed afterwards even on failure

### Scenario 6: Layout and project persistence at peak configuration

The persistence tail runs at the peak configuration reached above — X=AGE,
Y=HEIGHT, bins 200, rectangle shape, inverted color scheme — with nothing
reverted before saving.

Steps:
1. Save the view layout.
2. Add a Scatter plot viewer to the same view (a foreign viewer).
3. Re-apply the saved layout.
4. Verify the view's viewer set equals the SAVED set: a Density Plot is
   present AND the later-added Scatter plot is absent.
5. Verify the restored Density Plot kept its configuration: X column AGE,
   Y column HEIGHT, bins 200, bin shape rectangle, Invert Color Scheme on.
6. Delete the probe layout.
7. Save the view as a project through the ribbon **Save** button, dismiss the
   Share dialog, Close All, then reopen the saved project.
8. Verify a Density Plot viewer is present in the reopened view AND its
   configuration (X/Y columns, bins, bin shape, inverted color scheme) equals
   the persisted values — a cross-session round-trip.
9. Delete the probe project.

Expected:
- Re-applying the saved layout restores the SAVED viewer set (Density Plot present, Scatter plot absent) with the configured x/y/bins/binShape/invertColorScheme
- Reopening the saved project restores a Density Plot with the same configuration
- The probe layout and project are deleted even when a verification fails — they never leak

## Automation notes

Setup: the viewer handle is
`const dp = grok.shell.tv.viewers.find(v => v.type === 'Density plot');` views
are closed via `grok.shell.closeAll()`. The viewer is added via the Toolbox
icon `[name="icon-density-plot"]` (a synthetic `.click()` works for the icon).

Scenario 1 — selector popup: the popup opens on a real (trusted CDP) click or
a synthetic `mousedown` on the `.d4-column-selector-column` label inside
`[name="div-column-combobox-x"]` / `[name="div-column-combobox-y"]` — a
synthetic `.click()` does NOT open it. Wait for
`.d4-column-selector-backdrop`, type the column name with real keyboard
events (the first key opens the search input), then ArrowDown + Enter
commits. The product read is `dp.props.xColumnName` / `yColumnName` plus the
selector's `.d4-column-selector-column` DOM label text. The GROK-16612
console guard must filter infra noise — the shared dev server continuously
logs unrelated `WebSocket ... 503` reconnect errors.

Scenarios 2-3 — pixel measurements: the plot has a single
`canvas[name="canvas"]` (bins, axes, color scale are all drawn there). Ink is
measured as the non-white pixel fraction of the canvas backing store,
settle-gated on `d4-viewer-rendered`. The reference numbers are recon
orientation values (bins=5 hexagon ~0.43 vs bins=200 ~0.08; invert with
rectangle bins 0.27 -> 0.81; linear vs logarithmic transform is modest —
prefer mean luminance over the non-white fraction) and MUST be re-calibrated
live at spec-authoring time. The rectangle/hexagon shape delta at a fixed bin
count is modest (~0.045 at bins=50 in recon), so the shape switch may
honestly fall back to the no-error floor.

Scenario 2 — bins=1 lower boundary: the on-viewer `<input type="range">` bin
slider clamps at `min=1`, so bins=1 is the true lower edge. At bins=1 the whole
plot is one filled cell, so the ink jumps well above the 50-bin grid: live recon
measured `settledPx` totals px1 ~300972 vs px50 ~158509 (delta ~142463, ~1.9x)
on the 537x935 backing store; the spec's 70000 delta floor keeps a ~2x safety
margin and MUST be re-calibrated if the viewer size changes. The revert to 50
is asserted directionally (px1 - px50back over the same floor), not as an exact
pixel restore.

Spec floor calibration (settledPx totals on the 537x935 backing store;
re-calibrate if the viewer size changes): bins 5 vs 200 — px5 ~256947 vs
px200 ~44113 (delta ~212834, ~5.8x), spec floor 100000 (~2x margin); Invert
Color Scheme — pxBefore ~46656 -> pxInvert ~490723 (delta ~444067), spec
floor 200000 (~2x margin).

Scenario 5 — degenerate zero-width bin range: the fixture is a numeric constant
column created with
`df.columns.addNewCalculated('DP_ZW_CONST_<ts>', '42')` (int, min == max == 42)
and removed with `df.columns.remove(...)` in a `finally` (the MatrixPlot fixture
pattern). A numeric constant column IS offered by the numeric-only selector, so
it is assigned to X through the same on-viewer selector flow as Scenario 1
(`pickColumnViaSelector`, no JS-API fallback). Live recon confirmed the
divide-by-zero guard holds: assigning the constant column and rendering
rectangle then hexagon threw nothing, the viewer root stayed attached, and each
shape painted a non-zero degenerate strip (rectangle ~7346 px, hexagon ~11367
px). The honest floor asserted is therefore: no new page/console errors, root
attached, and a settled non-zero `settledPx` for both shapes. Teardown restores
X=AGE and binShape=rectangle before dropping the fixture column, so the
persistence tail still runs at the intended peak configuration.

Property-panel drives: `[name="prop-bins"]` — focus the
`property-grid-slider-textbox`, real Ctrl+A, type the value, Enter;
`[name="prop-bin-shape"]` and `[name="prop-color-transform-type"]` — click
the row's `[name^="prop-view-"]` cell, set the spinner select's `.value`,
dispatch `change`, Enter; `[name="prop-invert-color-scheme"]
input[type="checkbox"]` — a synthetic `.click()` toggles and repaints.
Collapsed property categories keep their rows in DOM (visibility-hidden);
JS-driven editor interaction works on hidden rows, but real-user clicks need
the category expanded first (e.g. click `[name="prop-category-misc"]`).

Scenario 4: the column names inside the popup grid are canvas-rendered (not
DOM text), so the numeric-only check reads the no-match behaviour — the
search stays open and Enter commits nothing. A no-match popup may not close
on Esc; close it with a synthetic `mousedown` on `document.body` (the
outside-click listener is on document mousedown).

Scenario 6: the layout is saved and re-applied via the JS API
(`tv.saveLayout()` / `grok.dapi.layouts.save` / `tv.loadLayout`) — the View >
Layout menu path has no headless handles; the round-trip end-state is the
same. The project is saved through the real ribbon Save button
(`[name="button-Save"]`) via the `saveProjectViaUI` helper from
`helpers/projects.ts`, because only the UI Save captures the view layout; the
"Share <project>" dialog that pops up after a successful save is dismissed via
its CANCEL button, and the save emits the benign console pair whitelisted per
canon. The probe layout and project names carry a `Date.now()` suffix so
concurrent runs never collide, and both are deleted in `finally` teardowns
(`grok.dapi.layouts.delete` / `grok.dapi.projects.delete`) so they are removed
even when an assertion fails.

---
{
  "order": 12,
  "datasets": ["System:DemoFiles/demog.csv"]
}
