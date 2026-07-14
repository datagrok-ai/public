# Bar chart tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Bar chart

## Stack column

1. Open Context Panel > **Stack** > set Stack to SEX — bars should split into stacked segments
2. Set Stack to RACE — stacked segments should update to race categories
3. Enable **Value > Relative Values** — all bars should become the same width
4. Disable **Relative Values** — bars return to absolute widths
5. Set Stack to None — stacking removed

## Sorting

1. Set Split to RACE
2. Open Context Panel > **Value > Bar Sort Type** > set to **by value** — bars should reorder by height
3. Set **Bar Sort Order** to **asc** — shortest bar first
4. Set **Bar Sort Order** to **desc** — tallest bar first
5. Set **Bar Sort Type** to **by category** — bars reorder alphabetically
6. Set **Bar Sort Order** to **asc** — A-Z order
7. Set **Bar Sort Order** to **desc** — Z-A order

## Value axis type

1. Set Value to HEIGHT, Split to RACE
2. Open Context Panel > **Value > Axis Type** > set to **log** — value axis switches to logarithmic scale
3. Set **Axis Type** back to **linear** — axis returns to normal
4. Set custom **Min** and **Max** values — bars should clip to the specified range
5. Clear Min and Max — auto range restored

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

## Aggregation types

1. Set Split to RACE, Value to AGE
2. Set **Value Aggr Type** to **avg** — bars show average AGE per category
3. Switch to **min** — bars show minimum values
4. Switch to **max** — bars show maximum values
5. Switch to **sum** — bars show sum values
6. Switch to **count** — bars show row counts
7. Set Value to WEIGHT, Aggr Type back to **avg** — bars show average weight

## Date/time split column

1. Set Split to STARTED — bars should group by date
2. Open Context Panel > **Category > Split Map** > set to **Year** — bars group by year
3. Change Split Map to **Month** — bars group by month
4. Change to **Quarter** — bars group by quarter

## Legend

1. Set Stack to SEX — legend should appear showing stack categories
2. Open Context Panel > **Legend > Legend Visibility** > set to **Always** — legend is always visible
3. Set **Legend Position** to different positions (Left, Right, Top, Bottom) — legend moves
4. Set Legend Visibility to **Never** — legend hides
5. Remove Stack column — legend should disappear automatically

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

## Orientation

1. Open Context Panel > **Orientation** > set to **horizontal** — bars render horizontally
2. Set to **vertical** — bars render vertically (rotated 90 degrees)
3. Set to **auto** — orientation chosen automatically based on viewer aspect ratio

## Data panel (SPGI dataset)

Setup: Close all, open demog and SPGI 

1. Go to the demog table
1. Add Bar chart
2. Open Context Panel > **Data > Row Source** — switch between All, Selected, Filtered
3. On the Context Panel > Data > Table, switch table to SPGI — bar chart should rebind to new table
5. Set **Filter** to `${CAST Idea ID} < 636500` — chart should show only matching rows
6. Set **Color > Color Column** to Chemical Space Y
7. Save layout via JS API
8. Close viewer
9. Apply the saved layout — filter, color coding, and all settings should restore
10. Delete the saved layout

## Filter Panel interaction

1. Open the Filter Panel
2. Apply a categorical filter — verify the bar chart updates
3. Apply a numeric filter — verify collaborative filtering between the Filter Panel and bar chart
4. Remove all filters

## Scrolling with range slider

1. Set Split to a column with many categories (e.g., Primary Scaffold Name)
2. Use the vertical range slider to scroll through the categories — verify bars scroll smoothly
3. Set Value to CAST Idea ID — verify scrolling still works with the new value column

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/demog.csv,System:DemoFiles/SPGI.csv"]
}
