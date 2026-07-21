---
feature: barchart
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes:
  - viewers.bar-chart
realized_as:
  - bar-chart-spec.ts
related_bugs: []
expected_results:
  - anchor: "Color coding"
    expectation: >-
      Setting Color Column to HEIGHT recolors the bars (canvas per-color delta
      well above the settle noise), and toggling Invert Color Scheme reverses the
      gradient with a second measurable recolor; Color Column reads back empty
      after reset.
  - anchor: "Include nulls"
    expectation: >-
      With a split column carrying missing values (HEIGHT, 751 nulls), toggling
      Include Nulls off removes the "missing" bar and back on restores it — each
      toggle produces a large canvas repaint delta.
  - anchor: "Bar style"
    expectation: >-
      Setting Bar Border Line Width 0 → 2 and Max Bar Height 50 → 20 each repaint
      the chart (canvas delta above the settle floor); the remaining style props
      read back their assigned values.
  - anchor: "Labels"
    expectation: >-
      Switching Show Labels from inside to never removes the in-bar value labels,
      changing the canvas; the other label modes read back correctly.
  - anchor: "Controls visibility"
    expectation: >-
      Turning the show*Selector props off collapses the in-chart value/category/
      stack selector nodes (computed display none, count 0) and turning them on
      lays them out again (count > 0).
  - anchor: "Aggregation types"
    expectation: >-
      Changing Value Aggr Type avg → max re-scales every bar and rebinding the
      value column AGE → WEIGHT re-scales again — both produce a measurable
      canvas delta; the props read back correctly.
  - anchor: "Legend position"
    expectation: >-
      Setting a Stack column materializes the legend host (absent before, present
      with ≥ 2 items after) and clearing the Stack column removes it again — a DOM
      round-trip; Legend Position reads back each of Left/Right/Top/Bottom.
  - anchor: "Title and description"
    expectation: >-
      With Show Title on, the title text "Demographics" appears in the viewer's
      panel titlebar and the description "By race" appears in the viewer body
      (DOM text asserts).
  - anchor: "Show values instead of categories"
    expectation: >-
      Enabling Show Values Instead Of Categories repaints the category axis strip
      (canvas delta above the settle floor); the prop reverts cleanly.
  - anchor: "Data panel"
    expectation: >-
      Row Source cycles through Filtered/Selected/All, the chart rebinds to the
      SPGI table, the filter and color-column apply, and after save → close →
      reload the layout restores both the color column and the filter (honest
      layout round-trip).
---

# Bar chart tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Bar chart

## Color coding

Actuation note: the recolor is asserted through the **Color Column** = HEIGHT
and **Invert Color Scheme** canvas repaints. The **Color Aggr Type** min/max/med
cycle is left read-back-only — each aggregation does recolor the bars, but the
per-aggr canvas delta is not separately asserted; the two anchored repaints carry
the visual proof.

1. Set Split to RACE, Value to AGE
2. Open Context Panel > **Color > Color Column** > set to HEIGHT — bars should be color-coded by average height
3. Change **Color Aggr Type** to **min** — colors should change based on minimum height per bar
4. Change to **max**, then **med** — verify colors update each time
5. Toggle **Invert Color Scheme** — color gradient should reverse
6. Set Color Column to None — color coding removed

## Include nulls

1. Set Split to a column that has missing values (HEIGHT — 751 nulls in demog;
   the categorical columns have none)
2. Open Context Panel > **Value > Include Nulls** — enabled by default, a "missing" bar should appear
3. Disable **Include Nulls** — the missing values bar should disappear
4. Re-enable — bar reappears

## Bar style

Actuation note (recon 2026-07-21): exact bar width/thickness is not
pixel-measurable headless, so the **Max Bar Height** and **Bar Border Line
Width** changes are asserted through the per-color canvas repaint delta
(`snapshotCanvasColors` + `diffCanvasColors`, with a settle-precheck to drain
the late render before measuring) rather than a geometric width read — the
delta is the drivable proxy for "the bars visibly changed". A precise width
check stays a human-side variation. **Bar Corner Radius**, **Vertical Align**
(Top/Bottom/Center), and **Show Category Zero Baseline** off are left
read-back-only: their on-canvas effect is a human-side variation, so the spec
covers them by property read-back rather than a canvas delta.

1. Open Context Panel > **Style**
2. Set **Bar Border Line Width** to 2 — bars should get visible borders
3. Set **Bar Corner Radius** to 10 — bar corners should become rounded
4. Set **Max Bar Height** to 20 — bars should become thinner
5. Change **Vertical Align** to Top, then Bottom, then Center — bars should align accordingly
6. Set **Show Category Zero Baseline** off — baseline line disappears

## Labels

Actuation note: the **inside → never** transition is asserted through the canvas
repaint (the in-bar label glyphs vanish). The **outside** and **auto** modes are
left read-back-only — the never transition carries the visual proof.

1. Set Split to RACE, Value to AGE
2. Open Context Panel > **Style**
3. Set **Show Labels** to **inside** — value labels appear inside bars
4. Set to **outside** — labels appear outside bars
5. Set to **never** — labels disappear
6. Set to **auto** — labels appear where they fit

## Controls visibility

Actuation note (recon 2026-07-21): the value/category/stack in-chart selectors
are DOM nodes whose computed `display` toggles (flex ↔ none) with the
show*Selector props — asserted by counting the laid-out selector nodes (exactly
3 when on, 0 when off; probed live). **Show Value Axis** and **Show Category
Values** are canvas-drawn — a live probe confirmed they toggle no DOM node (0
axis/category nodes) and produce a ~13 px canvas delta ≈ settle noise, so they
stay covered by the property read-back; their on-canvas appearance is a
human-side variation.

1. Open Context Panel
2. Toggle **Show Value Selector** off — value dropdown disappears from viewer
3. Toggle **Show Category Selector** off — category dropdown disappears
4. Toggle **Show Stack Selector** off — stack dropdown disappears
5. Toggle **Show Value Axis** off — value axis disappears
6. Toggle **Show Category Values** off — category labels disappear
7. Re-enable all — everything reappears

## Aggregation types

1. Set Split to RACE, Value to AGE
2. Set **Value Aggr Type** to **max** — bars show maximum AGE per category
3. Set Value to WEIGHT — bars rebind to the new value column

## Legend position

1. Set Stack to SEX and Legend Visibility to **Always** — legend appears
2. Set **Legend Position** to Left, Right, Top, Bottom — legend moves to each position
3. Remove Stack column — legend disappears automatically

## Title and description

1. Open Context Panel > **Description > Show Title** — enable
2. Set **Title** to "Demographics" — title appears at the top
3. Set **Description** to "By race" — description text appears
4. Change **Description Position** (Top, Bottom, Left, Right) — description moves
5. Set **Description Visibility Mode** to **Never** — description hides

## Show values instead of categories

1. Set Split to RACE, Value to AGE, Aggr Type to **avg**
2. Open Context Panel > **Category > Show Values Instead Of Categories** — category labels switch to show aggregated values
3. Disable — category names return

## Data panel

Setup: Close all, open both demog and full SPGI (System:DemoFiles/SPGI.csv)

1. Go to the demog table
2. Add Bar chart
3. Open Context Panel > **Data > Row Source** — switch between All, Selected, Filtered
4. On the Context Panel > Data > Table, switch table to SPGI — bar chart should rebind to new table
5. Set **Filter** to `${CAST Idea ID} < 636500` — chart should show only matching rows
6. Set **Color > Color Column** to Chemical Space Y
7. Save layout via JS API
8. Close viewer
9. Apply the saved layout — filter, color coding, and all settings should restore
10. Delete the saved layout

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/demog.csv,System:DemoFiles/SPGI.csv"]
}
