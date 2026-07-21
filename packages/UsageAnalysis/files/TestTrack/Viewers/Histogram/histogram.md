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

## Labels

1. Set Split Column to SEX
2. Set **Legend Visibility** to Never
3. Set **Legend Visibility** to Always
4. Change **Legend Position** to Top Right
5. Set Split Column to None
6. Enable **Description > Show Title**
7. Set **Title** to "Age Distribution"
8. Set **Description** to "Shows distribution of patient ages"
9. Set **Description Visibility Mode** to Always
10. Change **Description Position** to Bottom

## Context menu

1. Right-click on the histogram area -- context menu should show "Show Filtered Out Rows" and "Selection" group
2. Expand **Selection** -- should show "Show Current Row", "Show Mouse Over Row", "Show Mouse Over Row Group"
3. Toggle "Show Filtered Out Rows" from the context menu
4. Right-click on the X axis area (enable X axis first) -- menu should show "Show X Axis", "Controls Font", "Axis Font"

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
2. Go to the Context Panel > Data, set Row Source to Selected — verify the histogram shows only selected rows
3. Set Row Source to All — verify the histogram restores to full data
4. Set **Filter** formula to `${AGE} > 40`
5. Clear filter formula
6. Open Context Panel > **Data > Table**
7. Switch table to demog -- histogram updates to demog table columns
8. Set Value to WEIGHT
9. Switch back to SPGI

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
