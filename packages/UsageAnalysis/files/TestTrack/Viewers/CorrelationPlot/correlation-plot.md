---
feature: correlationplot
realizes_atlas:
  - correlationplot.cp.property-surface-smoke
  - correlationplot.int.menu-toggles-mirror-props
realizes:
  - viewers.correlation-plot
priority: p2
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-16818
    status: fixed
  - id: GROK-9310
    status: fixed
  - id: GROK-19052
    status: fixed
  - id: GROK-18487
    status: fixed
  - id: GROK-12586
    status: fixed
precondition_guards: []
expected_results:
  - anchor: "Title, description, and back color"
    expectation: >-
      After enabling Show Title and setting a title, the title text appears in
      the viewer's DOM; the description text appears while Description
      Visibility Mode is Always and disappears when it is Never; changing Back
      Color repaints the viewer background (a margin-point pixel probe or
      settle-gated canvas diff confirms the color change) with no errors raised.
  - anchor: "Misc property sequence and clean console"
    expectation: >-
      Switching the misc and style properties in sequence (Show Tooltip off/on,
      Ignore Double Click on/off, Col Header Font and Default Cell Font changed
      and restored) completes with a clean console — no errors are raised during
      the whole sequence (GROK-16818); the benign 'Package GrokML is not
      available...' line is filtered through a window-gated allowlist without
      masking other console errors.
  - anchor: "Context-menu toggles mirror properties"
    expectation: >-
      Clicking the Show Pearson R context-menu item flips the Show Pearson R
      setting to off AND the matrix visibly repaints (settle-gated canvas diff —
      the in-cell digits disappear and the value columns narrow); clicking
      Tooltip > Visible flips the Show Tooltip setting to off AND hovering a
      correlation cell afterwards produces no tooltip; re-clicking both menu
      items restores the settings and the behavior (round-trip).
  - anchor: "Pinned row-header columns"
    expectation: >-
      With the matrix wider than the viewer, the pinned type/name columns stay
      at the left edge — hovering the same row-header x still yields the
      column-statistics tooltip (the hard pinning assert). The wheel-driven
      column shift under a fixed x (a probe click reporting a different column
      pair after the horizontal scroll) is asserted when the wheel drives the
      scroll and is otherwise a classified reduction (status: waived,
      waiver_class: gesture-uncontrollable-headless; evidence: live-MCP recon —
      the horizontal scroll is a canvas-drawn range selector with no DOM thumb
      and no JS-readable scroll offset).
  - anchor: "Order or Hide Columns dialog"
    expectation: >-
      The Grid > Order or Hide Columns... dialog opens; its column filter select
      is driven through visible, hidden, and all and the column search field
      accepts and clears a query (real dialog controls — each drive refilters
      the embedded column list); the dialog closes cleanly via its CLOSE button,
      all with no exception (GROK-9310 regression guard), and the pinned name
      column stays in place (a row-header hover still shows the
      column-statistics tooltip). The per-column checkbox toggle is waived
      (status: waived, waiver_class: canvas-webgl-render) — the column list is
      an embedded canvas grid with no per-column DOM checkboxes
      (probe-verified).
  - anchor: "Grid color-coding apply to text"
    expectation: >-
      Enabling 'Apply to text' under Grid > Grid Color Coding on a participating
      column completes without errors — an honest no-error floor: the colored
      text is canvas-drawn and exposes no readable channel (GROK-19052
      regression guard); disabling it again also completes without errors.
  - anchor: "Table switch"
    expectation: >-
      After switching the Table property to spgi-100, the matrix recomputes for
      the new table without errors: getCorrelation on a spgi-100 numeric pair
      equals the runtime-computed Pearson reference over the same row mask
      within tolerance 1e-3, and the console stays clean (GROK-18487 regression
      guard).
  - anchor: "NaN edge cell"
    expectation: >-
      With a constant calculated fixture column participating in the matrix,
      hovering its correlation cell shows a tooltip whose text reads 'N/A'
      instead of a number, and no 'Unsupported operation' error appears
      (GROK-12586 regression guard); the fixture column is removed in teardown
      even on failure.
  - anchor: "Color-coding probes"
    expectation: >-
      The center pixels of three cells whose runtime-computed r values are
      positive, negative, and near zero fall into different hue families, with
      the near-zero cell visibly lighter than both others — point-wise per-cell
      pixel probes, never a global non-white count.
  - anchor: "Diagonal histograms repaint"
    expectation: >-
      A settle-gated canvas diff of a diagonal (histogram) cell region exceeds
      the settle noise after a filter change, and again after a row selection is
      made and after it is cleared — the diagonal histogram is rebuilt with the
      combined filter and the selection.
  - anchor: "Column width drag"
    expectation: >-
      A real mouse drag of a column-header edge repaints the matrix — the
      settle-gated canvas diff exceeds the settle noise; dragging the edge back
      to the original width also produces a canvas diff exceeding the settle
      noise (round-trip repaint); if the drag proves inert headless, a
      documented per-item reduction is recorded instead (gridCol.width is not
      JS-readable, so the diff is the only honest signal).
realized_as:
  - correlation-plot-spec.ts
---

# Correlation plot — property surface smoke

## Purpose

Verifies the Correlation plot's remaining property surface: title and description
text, back color, the misc-property clean-console sequence, context-menu toggles
that mirror Look properties, pinned row-header columns, the Order or Hide Columns
dialog, grid color-coding 'apply to text', switching the bound table, the
NaN/undefined-coefficient edge, per-cell color-coding probes, diagonal histogram
repaints, and column width dragging.

Matrix values, correlation type, row source, and persistence belong to
`correlation-plot-matrix-values.md`; cell gestures (click, double-click, hover
tooltip, Open as Table, Tooltip Edit/Properties) belong to
`correlation-plot-cell-interactions.md` — this scenario does not repeat them.

All sections start with the following sequence:
1. Close all open tables and viewers.
2. Open the **demog** dataset.
3. Add a Correlation plot by clicking its icon in the Toolbox > Viewers section.

## Title, description, and back color

1. On the Context Panel > Description, enable **Show Title** and set **Title**
   to "Correlation Analysis" — the title text appears on the viewer.
2. Set **Description** to "Shows pairwise correlations" and **Description
   Visibility Mode** to Always — the description text appears.
3. Change **Description Position** to Bottom.
4. Set **Description Visibility Mode** to Never — the description disappears.
5. On the Context Panel > Style, set **Back Color** to a light gray — the
   viewer background repaints.

Expected:
- Title and description text appear and disappear in the viewer as driven
  (Steps 1, 2, 4).
- The Back Color change is confirmed by a margin-point pixel probe or a
  settle-gated canvas diff, with no errors raised (Step 5).

## Misc property sequence and clean console

1. On the Context Panel, toggle **Show Tooltip** off, then on.
2. Toggle **Ignore Double Click** on, then off.
3. On the Context Panel > Style, set **Col Header Font** to a larger bold
   size, then restore it.
4. Set **Default Cell Font** to a larger size, then restore it.
5. Verify the console collected no errors during the whole sequence.

Expected:
- The whole switching sequence raises no console errors (GROK-16818); the
  benign GrokML availability warning is allowlisted for this window only.

## Context-menu toggles mirror properties

1. Right-click a correlation cell and click the **Show Pearson R** menu item.
2. Verify the **Show Pearson R** setting is now off and the matrix repainted
   — the in-cell digits disappear and the value columns narrow (settle-gated
   canvas diff).
3. Right-click again and click **Tooltip** > **Visible**.
4. Verify the **Show Tooltip** setting is now off, and hovering a
   correlation cell no longer shows a tooltip.
5. Restore both settings through the same menu items: **Show Pearson R**,
   then **Tooltip** > **Visible**.
6. Verify **Show Pearson R** and **Show Tooltip** are on again and the hover
   tooltip is back.

Expected:
- Each menu toggle flips the matching Look property AND produces its
  structural effect — the menu is the actuation path, the setting plus the
  repaint or tooltip behavior is the signal (round-trip, Steps 2, 4, 6).

## Pinned row-header columns

1. Shrink the viewer width until the matrix is wider than the viewer and the
   horizontal scrollbar appears.
2. Probe-click a value cell at a fixed x and record the column pair the
   cell-click event reports.
3. Wheel-scroll the matrix horizontally.
4. Probe-click at the same fixed x — the event reports a different column
   (the value columns scrolled).
5. Hover the pinned name column at its original x — the column-statistics
   tooltip still appears (the pinned type/name columns did not move).
6. Restore the viewer width.

Expected:
- The value columns scroll while the pinned type/name columns stay at the
  left edge (Steps 4–5) — structural pinning confirmed by the tooltip
  remaining on the row-header.

## Order or Hide Columns dialog

1. Right-click a correlation cell and open **Grid** > **Order or Hide
   Columns...**.
2. In the dialog, switch the column filter select through visible, hidden,
   and all — the embedded column list refilters on each switch.
3. Type "AGE" into the dialog's search field, then clear it.
4. Close the dialog with its CLOSE button.
5. Verify no exception was thrown and the pinned name column is still in
   place (a row-header hover still shows the column-statistics tooltip).

Expected:
- The dialog opens, the filter select and the search field drive through the
  real dialog controls, the dialog closes cleanly, and no exception is
  thrown; the name column stays in place (GROK-9310). The per-column
  checkbox toggle is waived — the column list is an embedded canvas grid
  with no per-column DOM checkboxes (probe-verified).

## Grid color-coding apply to text

1. Right-click the viewer and open **Grid** > **Grid Color Coding**.
2. Enable **Apply to text**.
3. Verify the action completes without errors.
4. Disable **Apply to text** again.
5. Verify the disable also completes without errors.

Expected:
- Both the enable and disable of 'Apply to text' complete without errors —
  honest no-error floor; the colored text is canvas-drawn with no readable
  channel (GROK-19052).

## Table switch

1. Open the spgi-100 dataset alongside demog.
2. On the Context Panel > Data, switch the **Table** property to spgi-100.
3. Verify the matrix recomputed for the new table: a spgi-100 numeric pair's
   correlation value equals the runtime-computed Pearson reference within
   tolerance 1e-3, and no errors appeared.

Expected:
- The viewer survives the table switch and recomputes correct values for
  the new table with a clean console (GROK-18487).

## NaN edge cell

1. Add a calculated column with a constant formula (e.g. `0`) to demog via
   Add New Column — its correlation with any column is undefined.
2. Ensure the fixture column participates in the matrix (add it to the Y
   columns if it did not join automatically).
3. Hover the fixture column's correlation cell — the tooltip shows 'N/A'
   instead of a number.
4. Verify no 'Unsupported operation' error appeared.
5. Remove the fixture column (teardown — runs even if a previous step
   failed).

Expected:
- The undefined-coefficient cell tooltip reads 'N/A' and no 'Unsupported
  operation' error is raised (GROK-12586); the fixture column is removed in
  teardown.

## Color-coding probes

1. Compute correlation values at runtime and pick three cells whose r values
   are clearly positive, clearly negative, and near zero.
2. Read the center pixel of each of the three cells.
3. Compare the three pixels: the positive-r and negative-r cells fall into
   different hue families, and the near-zero cell is visibly lighter than
   both.

Expected:
- Three point-wise per-cell probes confirm the [-1, 1] color scale: distinct
  hue families for opposite signs, light color near zero — never a global
  non-white count.

## Diagonal histograms repaint

1. Capture a settle-gated canvas snapshot of a diagonal (histogram) cell
   region as the baseline.
2. Apply a filter that hides part of the rows (e.g. narrow AGE on the
   Filter Panel) — the diagonal histogram is rebuilt; the settle-gated diff
   against the baseline exceeds the settle noise.
3. Clear the filter and re-baseline.
4. Select a band of rows — the diagonal cell repaints (diff exceeds settle
   noise).
5. Clear the selection — the diagonal cell repaints again.

Expected:
- The diagonal histogram cell repaints when the filter changes and when the
  selection is made and cleared (settle-gated per-cell canvas diffs,
  Steps 2, 4, 5).

## Column width drag

1. Position the mouse on the right edge of a value column's header.
2. Perform a real mouse drag to widen the column.
3. Verify the matrix repainted — the settle-gated canvas diff against the
   pre-drag baseline exceeds the settle noise.
4. Drag the edge back to approximately the original width.
5. Verify the matrix repainted again after the column width was restored —
   the settle-gated canvas diff exceeds the settle noise.

Expected:
- Both the widening drag and the restoring drag repaint the matrix (settle-gated
  canvas diffs, Steps 3 and 5); if the drag proves inert headless, record a
  documented per-item reduction — gridCol.width is not JS-readable, so the
  diff is the only honest signal.

## Automation notes

- Tooltip read channel: the single reusable .d4-tooltip element hosts the
  correlation tooltip content inline — its textContent carries the '<Type> R:
  <value>' line (probe-verified). An EMPTY .d4-tooltip element means the hover
  missed the cell (the element's idle state), not that the text lives on a
  hidden layer.
- Geometry staleness: re-read the viewer root's bounding rect immediately
  before every hover/click group that follows a viewer add/remove — dock
  relayout shifts the root, so cached Setup-time coordinates miss the 40x20
  cells. After hovering, poll for the '<Type> R:' line in the .d4-tooltip
  textContent instead of using fixed waits.
- Emission discipline for expected_results_coverage: every realized_by MUST be
  a verbatim, contiguous quote copied from the emitted spec (a single line is
  enough — e.g. the exact expect(...) line such as "expect(after).toBe(before)"
  with its real variable names); never paraphrase, never join multiple
  statements with ';', never use '(...)' or '…' ellipses — the mechanical
  verbatim matcher rejects them and Gate E fails with E-EXPECT-COVERAGE-01.
  For multi-line assertions quote only the final expect line.
- GrokML console noise: 'Package GrokML is not available...' is benign —
  filter it through a WINDOW-GATED allowlist (allowlist only within the
  observed sequence, never a blanket filter) WITHOUT masking other console
  errors (GROK-16818).
- Cell and menu coordinates: compute cell centers from the viewer root's
  bounding rect plus the documented grid arithmetic (pinnedW / headerH /
  cellW 40-or-20 / rowH 20 — refdoc §Cell Geometry and Mouse Clicks),
  calibrated by one probe click whose event reports the pair. Context-menu
  items have no name= attributes — find by .d4-menu-item-label text with the
  closest('.d4-menu-item') identity check (a group's own label comes AFTER
  its nested children in the popup DOM).
- Diagonal-selection diff actuation: the selection-made/cleared diff runs
  with Row Source switched to Selected for the measurement window (restored
  after) — a documented actuation deviation: under the default Filtered row
  source a bare selection does not repaint the diagonal, so asserting there
  would be a false failure. The body step's user-visible claim (the histogram
  repaints on selection) is unchanged.
- Prop mirror reads and diff channels: the context-menu toggle steps read
  the flip back via props.showPearsonR / props.showTooltip (the menu is the
  actuation path; the prop plus its structural effect is the signal); the
  Table switch step actuates via props.table (the Context Panel > Data
  Table selector writes the same property); the diagonal-histogram repaint
  is driven through df.filter and the row selection — these channels live
  in the spec, never in the step prose.
- Spec must keep (Gate E checklist): the restore-drag snap/diff assert
  (width-drag round-trip); the setToText driven-guard (the apply-to-text
  write proven, not just floored); the pinned-column tooltip matched by
  /avg:|min:/; descriptionPosition covered with a DOM-position assert;
  the visibility-based off-signal for the disabled tooltip (the reusable
  .d4-tooltip keeps stale text — computed display is the honest channel);
  explicit null returns from evaluates touching DG objects (implicit
  Dart-backed returns fail serialization); per-section recalibration after
  any view add/close; menu items resolved by own-label TEXT, never name=
  attributes (grid-generic items carry no name).
- Spec imports: this scenario lives in the nested Viewers/CorrelationPlot/
  section — spec imports must be '../../spec-login' and '../../helpers/*'.

---
{
  "order": 10,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
