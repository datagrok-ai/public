---
feature: histogram
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes:
  - viewers.histogram
realized_as:
  - histogram-spec.ts
related_bugs: []
expected_results:
  - anchor: "Spline mode"
    expectation: >-
      Driving spline and fill-spline on and off raises no console or page errors
      and leaves the histogram alive (canvas-only, no-error floor).
  - anchor: "Appearance — selector visibility"
    expectation: >-
      Turning Show Column Selector off hides the value-column combobox and turning
      it back on restores it; Show Split Selector does the same for the split-column
      combobox — each control's DOM visibility round-trips with its property.
  - anchor: "Appearance — canvas floor"
    expectation: >-
      Toggling the axis visibility props, X-axis height, Allow Column Selection,
      Show Bin Selector and Show Range Slider and reverting to defaults raises no
      console or page errors and leaves the histogram alive (canvas-only, no-error
      floor).
  - anchor: "Labels Step 1"
    expectation: >-
      Splitting by SEX renders a legend with exactly the two categories F and M;
      splitting by RACE renders four legend items — the legend item count tracks
      the split column's categories.
  - anchor: "Labels Step 2"
    expectation: >-
      Setting Legend Visibility to Never removes the legend host from the DOM;
      setting it back to Always restores it.
  - anchor: "Labels Step 6"
    expectation: >-
      Setting a Description makes its text appear inside the viewer's rendered
      content; clearing it removes the text.
  - anchor: "Context menu Step 1"
    expectation: >-
      Right-clicking the histogram canvas opens a menu that contains the
      histogram entries Show Filtered Out Rows, Selection, Show Current Row, Show
      Mouse Over Row, Show Mouse Over Row Group, plus the axis entries Show X
      Axis, Axis Font, Controls Font.
  - anchor: "Context menu Step 3"
    expectation: >-
      Clicking Show Filtered Out Rows in the context menu flips the viewer's
      showFilteredOutRows property (menu-to-prop round-trip).
  - anchor: "Layout persistence Step 8"
    expectation: >-
      After saving the layout, closing the viewer, and re-applying the layout,
      the histogram restores with WEIGHT value, 15 bins, RACE split, and stacked
      mode.
  - anchor: "Data properties Step 2"
    expectation: >-
      Setting Row Source to Selected shrinks the histogram's rendered canvas
      content against the All baseline and removes the orange selection overlay
      (present under All, zero under Selected); switching back to All restores both
      the content and the overlay exactly.
  - anchor: "Data properties Step 4"
    expectation: >-
      Setting the Filter formula to ${AGE} > 40 shrinks the histogram's rendered
      canvas content; clearing the formula restores it exactly.
  - anchor: "Data properties Step 7"
    expectation: >-
      Switching the viewer's Table from SPGI to demog rebinds the histogram to
      the demog dataframe and re-picks a value column that belongs to demog.
---

# Histogram tests (Playwright)

> **Note:** Bins configuration, Bin selection, Split column, Color coding, Value range, and
> Filtering are now covered by the focused scenarios `histogram-setup-tune-select.md`,
> `histogram-range-filter-bounds.md`, and `histogram-split-and-color.md`. This scenario retains
> only the areas without a focused cp yet (Spline mode, Appearance, Labels, Context menu,
> Layout persistence, Data properties).

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Histogram

## Spline mode

1. Enable **Spline**
2. Enable **Fill Spline**
3. Disable **Fill Spline**
4. Disable **Spline**

Actuation note (recon 2026-07-21): spline and fill-spline are pure canvas
painting — there is no DOM element or dataframe value that changes when they
toggle. The automated check is therefore the no-error floor (driving all four
toggles raises no console/page error and the viewer stays alive), not a
property read-back. A visual before/after comparison stays a human variation.

## Appearance

1. Enable **Value > Show X Axis**
2. Enable **Misc > Show Y Axis**
3. Disable both
4. Set **Style > X Axis Height** to 30
5. Disable **Show Column Selector**
6. Disable **Show Bin Selector**
7. Disable **Show Split Selector**
8. Disable **Show Range Slider**
9. Re-enable all four

Actuation note (recon 2026-07-21): two of these toggles DO carry a real DOM
signal. **Show Column Selector** drives the visibility of the value-column
combobox (`[name="div-column-combobox-value"]`) and **Show Split Selector**
drives the split-column combobox (`[name="div-column-combobox-split"]`, visible
only while a Split Column is set). A live probe flipped each prop off→on three
times and the combobox computed-visibility tracked it deterministically (off →
hidden, on → visible). Those two are asserted as a real present→absent→present
round-trip. **Show Bin Selector**, **Show Range Slider**, the axis toggles and
X-axis height have no headless DOM signal — the axes and the bin/range strip are
canvas-drawn — so they ride the no-error floor (drive every toggle, revert to
defaults, no error, viewer alive). Verifying the bin/range controls actually
appear/hide is a human variation.

## Labels

1. Set Split Column to SEX, then RACE -- the legend lists exactly the split
   column's categories (SEX → F, M; RACE → four items)
2. Set **Legend Visibility** to Never (legend host removed), then Always
   (restored)
3. Change **Legend Position** to Top Right (legend still rendered)
4. Set Split Column to None
5. Enable **Description > Show Title**
6. Set **Title** to "Age Distribution", **Description** to "Shows distribution
   of patient ages", **Description Visibility Mode** to Always, and **Description
   Position** to Bottom -- the description text appears inside the rendered
   viewer
7. Clear title and description -- the description text is gone

Actuation note (recon 2026-07-21): the legend is a real DOM signal
(`[name="legend"] .d4-legend-item` count tracks the split categories) and the
description is read back from the viewer's rendered `innerText`. The **title**
lives in the surrounding header chrome (not inside the viewer element), so it is
driven but not asserted — same treatment as PC Plot. Legend **position** is a
canvas/layout property with no per-position DOM signal, so it is driven only to
confirm the legend survives the layout pass.

## Context menu

1. Right-click on the histogram canvas -- the context menu shows the histogram
   entries "Show Filtered Out Rows" and the "Selection" group, plus the axis
   entries "Show X Axis", "Axis Font", "Controls Font"
2. The "Selection" group offers "Show Current Row", "Show Mouse Over Row", "Show
   Mouse Over Row Group"
3. Toggle "Show Filtered Out Rows" from the context menu -- the viewer's
   `showFilteredOutRows` property flips
4. (Human variation) Right-click specifically on the X axis area

Actuation note (recon 2026-07-21): right-clicking the histogram canvas opens a
single flattened menu (`.d4-menu-item-label`) that already contains BOTH the
table/view entries AND the histogram-specific entries — the old spec's claim
that only the view menu appears was wrong. The automated check asserts the
promised histogram items are present and that clicking "Show Filtered Out Rows"
round-trips the property. Step 4's position-precise right-click on the X-axis
region is reduced to the same axis items ("Show X Axis", "Axis Font", "Controls
Font") that the main-area menu already surfaces, because the axis strip is
canvas-drawn with no DOM handle to target headless; the axis-region right-click
stays a human variation.

## Layout persistence

1. Set Value to WEIGHT
2. Set Bins to 15
3. Set Split Column to RACE
4. Enable **Split Stack**
5. Save the layout via JS API
6. Close histogram viewer
7. Apply the saved layout
8. Verify histogram restores with WEIGHT, 15 bins, RACE split, and stacked mode

## Data properties

Setup: open SPGI dataset, add Histogram.

1. Select some rows in the grid
2. Go to the Context Panel > Data, set Row Source to Selected — the histogram's
   rendered canvas content shrinks against the All baseline
3. Set Row Source to All — the rendered content restores to the full-data
   baseline exactly
4. Set **Filter** formula to `${AGE} > 40` — the rendered canvas content shrinks
5. Clear filter formula — the rendered content restores exactly
6. Open Context Panel > **Data > Table**
7. Switch the viewer's Table to demog — the histogram rebinds to the demog
   dataframe and re-picks a demog value column
8. Switch back to SPGI

Actuation note (recon 2026-07-21):

- **Row Source** (Step 2-3) is measured two independent ways. The non-white
  canvas total shrinks when Selected paints only the selected rows (All 433634 px
  → Selected 430846 px, Δ 2788), restoring exactly when switched back (restore
  noise 0 across selection sizes 3–500). The orange selection-overlay hue
  (`countSelectionHuePixels`) is the second signal: 566 px under All (the overlay
  highlights the selected rows on top of the full-data bars), exactly 0 under
  Selected (no overlay when the whole chart IS the selection), 566 px restored.
  SPGI is used because it carries a selection large enough to measure. The prior
  spec only read `props.rowSource` back (a vacuous set-then-read), which this
  replaces.
- **Filter formula** (Step 4-5) is measured the same way (cleared 229645 px →
  filtered 228276 px, restored exactly). The viewer `filter` prop filters only
  which rows the histogram paints — it does NOT change
  `df.filter.trueCount` — so the dataframe filter count is not a valid signal
  here. **demog** is used because the filter expression references its **AGE**
  column, which SPGI does not have; on SPGI the same formula would be inert.
- **Table switching** (Step 7) is implemented via `props.table = 'demog'` with
  both tables open, asserting the bound `dataFrame.name` changes SPGI → demog and
  the auto-picked value column becomes a demog column (CAST Idea ID → AGE). This
  replaces the prior silent SKIP.

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
