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
---

# Bar chart tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Bar chart

## Color coding

1. Set Split to RACE, Value to AGE
2. Open Context Panel > **Color > Color Column** > set to HEIGHT — bars should be color-coded by average height
3. Change **Color Aggr Type** to **min** — colors should change based on minimum height per bar
4. Change to **max**, then **med** — verify colors update each time
5. Toggle **Invert Color Scheme** — color gradient should reverse
6. Set Color Column to None — color coding removed

## Include nulls

1. Set Split to a column that has missing values
2. Open Context Panel > **Value > Include Nulls** — enabled by default, a "missing" bar should appear
3. Disable **Include Nulls** — the missing values bar should disappear
4. Re-enable — bar reappears

## Bar style

1. Open Context Panel > **Style**
2. Set **Bar Border Line Width** to 2 — bars should get visible borders
3. Set **Bar Corner Radius** to 10 — bar corners should become rounded
4. Set **Max Bar Height** to 20 — bars should become thinner
5. Change **Vertical Align** to Top, then Bottom, then Center — bars should align accordingly
6. Set **Show Category Zero Baseline** off — baseline line disappears

## Labels

1. Set Split to RACE, Value to AGE
2. Open Context Panel > **Style**
3. Set **Show Labels** to **inside** — value labels appear inside bars
4. Set to **outside** — labels appear outside bars
5. Set to **never** — labels disappear
6. Set to **auto** — labels appear where they fit

## Controls visibility

1. Open Context Panel
2. Toggle **Show Value Selector** off — value dropdown disappears from viewer
3. Toggle **Show Category Selector** off — category dropdown disappears
4. Toggle **Show Stack Selector** off — stack dropdown disappears
5. Toggle **Show Value Axis** off — value axis disappears
6. Toggle **Show Category Values** off — category labels disappear
7. Re-enable all — everything reappears

## Aggregation types (max + value-column switch)

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

## Data panel (demog → SPGI table switch)

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
