# Scatter Plot — Test Coverage Explanation

> This document explains **why** each section in `Scatter plot - additional tests.md` is not covered
> by the existing test artifacts and what gap each section fills.

---

## Existing Coverage Baseline

### Manual test track (`Scatter plot.md`)
Covers: table switching, zoom/filter modes, axis min/max, log scale, axis inversion, coloring (linear/categorical/conditional), size-coding, marker types (categorical/structure/datetime columns), jitter, selection (Shift+drag, Ctrl+click, lasso), tooltip, formula lines, structures as axes, filter panel interaction, empty column on log scale, negative/zero values on log scale, legend color picker.

### Manual test track (`Legend\Scatterplot.md`)
Covers: combined color+marker legend, legend update on axis change, in-viewer filtering via expression, filter panel → legend sync, color coding from grid column.

### Automated tests (`scatter-plot.ts`)
Covers four regression bug fixes only:
- Wrong axis range (#1613, #1858)
- Log scale with near-zero values (#1764)
- Table unexpectedly filtered after empty column creation (#1744)

### Automated tests (`rendering-scatter-plot-tooltip-tests.ts`)
Covers: molecule tooltip rendering on hover via `mousemove` event simulation.

### Automated tests (`scatter-plot-3d.ts`)
Covers: 3D scatter plot only (visual, API, serialization).

---

## Gap Analysis by Section

---

### Section 1 — Axis Histograms

**Why not covered:**
`scatterplot_look.dart` exposes `showXHistogram`, `showYHistogram`, and `histogramBins` as first-class properties, but no manual test step and no automated test references these properties. The existing `Scatter plot.md` does not mention histograms at all. This is a distinct visual sub-component rendered alongside the main chart and deserves its own verification pass.

**What the tests add:**
- Verify the histogram appears and disappears when toggled.
- Verify `histogramBins` changes the number of bars.
- Verify the histogram responds to active row filters (not just static data).
- Negative: confirm no crash or silent failure when the feature is toggled on a categorical axis (where it has no defined behavior).

---

### Section 2 — Error Bars / Whiskers

**Why not covered:**
The `_ScatterPlotWhiskersFeature` class and the six whisker properties (`xWhiskerMinColumnName`, `xWhiskerMaxColumnName`, `xWhiskerRangeColumnName`, and Y equivalents) are entirely absent from both manual tracks and all automated tests. Auto-detection logic (`findErrorColumnsForAxis`) and the mutual-exclusion rule between Range and Min/Max are non-trivial behaviors with no test coverage whatsoever.

**What the tests add:**
- Verify auto-detection by column naming convention.
- Verify manual assignment of Min/Max and Range columns.
- Verify the mutual-exclusion contract (setting Min clears Range).
- Negative: confirm that whiskers with negative values on a log scale do not produce rendering artifacts or exceptions.

---

### Section 3 — Lines Order

**Why not covered:**
`_ScatterPlotLinesOrderFeature` and the `linesOrderColumnName` / `linesWidth` properties are not mentioned anywhere in the manual test track or in any automated test. This feature changes the chart appearance significantly (connecting markers with lines per category) and its interaction with the `colorCol` path is a distinct rendering code path not exercised elsewhere.

**What the tests add:**
- Verify line drawing order follows column values.
- Verify per-category line coloring when a categorical Color column is set.
- Verify `linesWidth` property has a visible effect.
- Negative: confirm graceful rejection of a non-numeric column for Lines Order.

---

### Section 4 — Regression Line Extended Checks

**Why not covered:**
The existing manual track only mentions enabling the regression line and saving to layout. It does not cover:
- The five optional statistics (Spearman, Pearson, MAE, RMSE, equation) individually.
- The tooltip/click-through behavior on each statistic label (Wikipedia links).
- The `regressionPerCategory` mode.
- The keyboard shortcut `R` for toggling.
- Regression behavior specifically on a logarithmic scale (a separate code path in `currentParams`).
- The degenerate case of a single visible data point.

**What the tests add:**
All of the above are added as discrete, individually verifiable steps, making regression line testing more exhaustive and reproducible.

---

### Section 5 — Smart Labels Extended Checks

**Why not covered:**
`_ScatterPlotSmartLabelsFeature` is the most complex feature module (~300+ lines) but receives no dedicated test coverage. The manual track mentions structures as axes (which triggers molecule label rendering) but does not address:
- Multiple label columns simultaneously.
- The four `showLabelsFor` modes (`MouseOverRow`, `MouseOverGroup`, `Selected`, `All`).
- The three `displayLabels` modes (`Always`, `Auto`, `Never`).
- `useLabelAsMarker` and `labelAsMarkerSize`.
- User-initiated label dragging and persistence.
- `showLabelNamedColumns` visibility setting.
- The guard that disables label dragging when the Dot Renderer is active (a code-level invariant: `grabbingAllowed = !dotRendererApplied && …`).

**What the tests add:**
Each of the above modes and interactions is given an explicit test step, covering both the happy path and the known limitation (drag disabled on Dot Renderer).

---

### Section 6 — Drop Lines

**Why not covered:**
`showDropLines` is defined in `scatterplot_look.dart` and implemented in `features/drop_line.dart`. Neither the manual track nor any automated test references this property. The feature is interactive (depends on mouse position) and has coordinate-system-specific behavior on inverted and logarithmic axes.

**What the tests add:**
- Basic enable/disable verification.
- Axis inversion interaction (coordinate displayed must match the inverted scale).
- Log scale interaction (coordinate must use the log-transformed value).

---

### Section 7 — Keyboard Navigation

**Why not covered:**
`_ScatterPlotNavigationFeature` registers handlers for arrow keys, `+`/`-`, `H`, `L`, `R`, and `M`. The manual track only mentions lasso (`Shift+L`) in passing. No test verifies arrow-key row navigation, zoom keys, or the `H` (Home/reset) shortcut. These are user-discoverable interactions that frequently regress silently.

**What the tests add:**
Explicit verification of the four most important keyboard interactions: row navigation, zoom, viewport reset, and lasso toggle. Note: `R` is also covered in Section 4.3 as part of the regression line context.

---

### Section 8 — Markers Additional Types and Edge Values

**Why not covered:**
The existing manual track sets markers to categorical, structure, and datetime columns (which changes which *column* drives marker shape) but does not test:
- Individual marker *type* values (circle, square, triangle, etc.) when `markersColumnName` is null.
- `markerDrawBorder` and `markerBorderWidth`.
- The edge case `markerMinSize == markerMaxSize` with size-coding active (potential division-by-zero).
- `markerOpacity = 0` (transparency with functional hit-testing).
- `markerMinSize > markerMaxSize` constraint.

**What the tests add:**
Each of these properties is verified in isolation, exposing potential rendering bugs and constraint violations.

---

### Section 9 — filterOutInvalid and showFilteredOutPoints

**Why not covered:**
`filterOutInvalid` and `showFilteredOutPoints` are separate Data-category properties. The existing manual track and bug-regression tests partially touch the null/empty column scenario (#1744), but they do not test:
- `filterOutInvalid = true` as a deliberate user setting (separate from the bug scenario).
- `showFilteredOutPoints = true` and its visual effect.
- The interaction between `showFilteredOutPoints` and `Zoom And Filter = filter by zoom` (do filtered-out points remain visible or get hidden?).
- Re-verification of #1744 (all-null column + filterOutInvalid) as a confirmed fix validation.

**What the tests add:**
Explicit positive verification of both properties and their interaction, plus a regression guard for #1744.

---

### Section 10 — Color Scale Extended Checks

**Why not covered:**
The manual track covers color-coding column assignment and the color picker in the legend. It does not cover:
- `colorAxisType = logarithmic` (a separate rendering path for the color gradient scale).
- `colorMin` / `colorMax` manual clamping (value range clipping on the color scale).
- Click-drag on the color scale bar to filter rows (an interactive gesture distinct from the legend).
- `invertColorScheme` (the visual inversion of the gradient).

**What the tests add:**
These four properties represent distinct rendering behaviors and user interactions on the color scale component, none of which are tested anywhere.

---

### Section 11 — Show Current / MouseOver / Selected Points

**Why not covered:**
`showCurrentPoint`, `showMouseOverPoint`, `showMouseOverRowGroup`, and `showSelectedRows` are four independent boolean properties in `scatterplot_look.dart`. The manual track verifies that selection works but never explicitly disables these visual indicators to verify the off-state. Disabling them affects multiple rendering layers independently and can mask bugs where the indicator is incorrectly suppressed or still drawn.

**What the tests add:**
Each property is toggled off and the resulting visual state is verified explicitly. This is a "control surface" test that guards against property changes having no effect.

---

### Section 12 — axesFollowFilter

**Why not covered:**
`axesFollowFilter` (default `true`) causes the viewport to auto-zoom when a range filter is applied from the Filter Panel. This is tested indirectly in the Zoom And Filter manual steps, but not for the specific case of `axesFollowFilter = false`, which is a user-visible setting in the Data properties panel. The off-state (viewport frozen despite filter) is never verified.

**What the tests add:**
Explicit verification that the default on-state zooms correctly, and that the off-state correctly freezes the viewport independent of the filter state.

---

### Section 13 — Axis Label Orientation

**Why not covered:**
`xAxisLabelOrientation` has four choices (`Auto`, `Horz`, `Vert`, `45 degrees`) declared in `scatterplot_look.dart`. Neither the manual track nor any automated test changes this property. It is particularly important for categorical axes with long names, where incorrect orientation produces overlapping or clipped labels.

**What the tests add:**
Systematic verification of all four orientations, plus an edge case with long category names to confirm the auto-fallback or truncation behavior.

---

### Section 14 — API / Serialization

**Why not covered:**
The automated API tests (`scatter-plot.ts`) only test four bug regressions via `testViewer`. They never verify:
- Multi-option `scatter()` calls with Color, Size, markerType, showRegressionLine, and zoomAndFilter set together.
- The `onAfterDrawScene` event and the `worldToScreen()` coordinate transform (covered in Chem only for molecule tooltips, not general coordinates).
- The `onZoomed` event (listed in `scatter-plot-events.js` sample but never asserted in tests).
- Layout save/restore for whisker and lines-order properties specifically.
- `setOptions()` with an invalid column name (error handling path).

**What the tests add:**
API surface coverage for events, coordinate transforms, and error handling — areas exercised in examples but not in assertions.

---

### Section 15 — Trellis Plot

**Why not covered:**
`scatterplot_meta.dart` declares `trellisProperties` (`xColumnName`, `yColumnName`, `colorColumnName`, `sizeColumnName`), confirming the scatter plot is a supported inner viewer for Trellis. Neither the manual track nor any automated test exercises the scatter plot inside a Trellis container. Trellis adds its own rendering loop and selection-sync mechanism that are orthogonal to the standalone viewer tests.

**What the tests add:**
Basic smoke test for scatter plot as a Trellis inner viewer, and verification that cross-panel selection synchronization works correctly.

---

### Section 16 — Performance / WebGPU Fallback

**Why not covered:**
`PowerGrid/src/webgpu/scatterplot.ts` implements the WebGPU renderer and `PowerGrid/src/tests/webgpu-tests.ts` benchmarks it at 10K / 1M / 10M rows, but these tests are marked as **benchmark-only** and skipped in normal test runs. No functional test verifies:
- That the Dot Renderer / WebGPU path activates automatically at the correct row count threshold (>50K).
- That `hitTest` (point click → row index) works correctly in the fast renderer.
- That the system falls back gracefully to Canvas rendering when WebGPU is unavailable.

**What the tests add:**
Functional (not benchmark) verification of renderer selection and fallback behavior under realistic conditions.

---

### Section 17 — Negative / Edge Case Scenarios

**Why not covered:**
The existing test track is entirely positive-path oriented. The automated tests cover four specific known-bug regressions. No test deliberately exercises:

| Scenario | Why it matters |
|----------|---------------|
| X Min > X Max | Constraint inversion — can produce an infinite rendering loop or inverted axis without visible error |
| Log scale on all-zero column | Division by zero / `log(0) = -Infinity` — can produce NaN in canvas coordinates |
| Very long label strings | Canvas `fillText` with unbounded strings can overflow the clip region |
| All regression statistics enabled | Formula box layout overflow — each statistic takes vertical space; all five together may exceed the viewer bounds |
| Shift+Drag in Pan mode | The Shift modifier overrides Mouse Drag mode; if not handled, selection becomes impossible in Pan mode |
| Empty dataframe (0 rows) | Null-guard paths throughout the rendering pipeline — a common source of crashes in viewer initialization |
| Delete column assigned to axis | Runtime column removal — the viewer holds a direct `Column` reference; removal must be handled gracefully |

**What the tests add:**
Systematic negative coverage that confirms the viewer fails safely and predictably under invalid or extreme inputs, rather than crashing or producing silent incorrect output.

---

## Summary Table

| Section | Feature / Property | Source file | Gap type |
|---------|--------------------|-------------|----------|
| 1 | showXHistogram, showYHistogram, histogramBins | scatterplot_look.dart | Not mentioned in any existing test |
| 2 | xWhiskerMin/Max/Range, yWhiskerMin/Max/Range | scatterplot_look.dart, features/whiskers.dart | Not mentioned in any existing test |
| 3 | linesOrderColumnName, linesWidth | scatterplot_look.dart, features/lines_order.dart | Not mentioned in any existing test |
| 4 | showSpearmanCorrelation, showPearsonCorrelation, showMAE, showRMSE, regressionPerCategory, R key | features/regression_line.dart | Partially mentioned; statistics, key shortcut, and log scale not tested |
| 5 | labelColumnNames, showLabelsFor, displayLabels, useLabelAsMarker, label drag | features/smart_labels.dart | Molecule labels mentioned; all other modes untested |
| 6 | showDropLines | scatterplot_look.dart, features/drop_line.dart | Not mentioned in any existing test |
| 7 | Arrow keys, +/-, H key | features/navigation.dart | Not mentioned in any existing test |
| 8 | markerType, markerDrawBorder, markerBorderWidth, markerOpacity, markerMinSize/MaxSize | scatterplot_look.dart | Column-driven markers tested; direct type/style properties not tested |
| 9 | filterOutInvalid, showFilteredOutPoints | scatterplot_look.dart | Only bug #1744 scenario partially covered |
| 10 | colorAxisType, colorMin, colorMax, invertColorScheme, color scale bar interaction | scatterplot_look.dart | Column assignment covered; scale properties not tested |
| 11 | showCurrentPoint, showMouseOverPoint, showMouseOverRowGroup, showSelectedRows | scatterplot_look.dart | Visual state when disabled never verified |
| 12 | axesFollowFilter | scatterplot_look.dart | Off-state never verified |
| 13 | xAxisLabelOrientation | scatterplot_look.dart | Not mentioned in any existing test |
| 14 | onAfterDrawScene, onZoomed, worldToScreen, setOptions error path | API surface | Events and error path not asserted |
| 15 | Trellis inner viewer, cross-panel selection | scatterplot_meta.dart | Not mentioned in any existing test |
| 16 | Dot Renderer / WebGPU activation, fallback | PowerGrid/src/webgpu/scatterplot.ts | Benchmarks exist but skipped; functional check missing |
| 17 | Constraint violations, empty data, column deletion | Multiple | No negative tests exist in the current track |

---

{
  "order": 11,
  "datasets": []
}
