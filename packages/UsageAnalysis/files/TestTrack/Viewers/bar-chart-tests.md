# Bar chart tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Bar chart

## Stack column

1. Open Context Panel > **Stack** > set Stack to SEX — bars should split into stacked segments by sex
2. Set Stack to RACE — stacked segments should update to race categories
3. Enable **Value > Relative Values** — all bars should become the same width, showing relative proportions
4. Disable **Relative Values** — bars return to absolute widths
5. Set Stack to None — stacking removed

## Sorting

1. Set Split to RACE
2. Open Context Panel > **Value > Bar Sort Type** > set to **by value** — bars should reorder by their height
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

## On click actions

1. Open Context Panel > **Data > On Click** > set to **Filter**
2. Click a bar — only the rows in that category should remain visible in the grid
3. Click empty space near the chart — filter should reset
4. Set **On Click** to **Select** — clicking a bar should select those rows in the grid
5. Set **On Click** to **None** — clicking should have no effect

## Include nulls

1. Set Split to a column that has missing values
2. Open Context Panel > **Value > Include Nulls** — enabled by default, a "missing" bar should appear
3. Disable **Include Nulls** — the missing values bar should disappear
4. Re-enable — bar reappears

## Show filtered rows

1. Open Context Panel > **Data > Row Source** > set to **All**
2. Open Filter Panel and filter some rows
3. Enable **Show Filtered Rows** — filtered portion should be indicated on bars
4. Disable **Show Filtered Rows** — indication disappears

## Bar style

1. Open Context Panel > **Style**
2. Set **Bar Border Line Width** to 2 — bars should get visible borders
3. Set **Bar Corner Radius** to 10 — bar corners should become rounded
4. Set **Max Bar Height** to 20 — bars should become thinner
5. Change **Vertical Align** to Top, then Bottom — bars should align accordingly
6. Set **Show Category Zero Baseline** off — baseline line disappears

## Labels

1. Set Split to RACE, Value to AGE
2. Open Context Panel > **Style** (or right-click > Properties)
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

## Orientation

1. Open Context Panel > **Orientation** > set to **horizontal** — bars render horizontally (default-like)
2. Set to **vertical** — bars render vertically (rotated 90 degrees)
3. Set to **auto** — orientation chosen automatically based on viewer aspect ratio
4. Resize the viewer to be very tall and narrow — in auto mode, bars should switch to vertical

## Aggregation types

1. Set Split to RACE, Value to AGE
2. Set **Value Aggr Type** to **avg** — bars show average AGE per category
3. Switch to **min** — bars show minimum values
4. Switch to **max** — bars show maximum values
5. Switch to **sum** — bars show sum values
6. Switch to **count** — bars show row counts (Value column selector should hide)
7. Set Value to WEIGHT, Aggr Type back to **avg** — bars show average weight

## Zoom and pan

1. **Alt + drag** vertically on the chart — chart should zoom into that range of categories
2. Use **mouse wheel** to scroll through categories when zoomed in
3. **Drag** on empty space (no modifiers) — chart should pan vertically through categories
4. **Double-click** the chart — view should reset to show all categories

## Selection interactions

1. Set **On Click** to **Select**
2. Click a bar — rows in that category are selected (highlighted in grid)
3. **Ctrl + click** another bar — selection adds that category
4. **Shift + drag** a rectangle over multiple bars — all covered categories are selected
5. Click empty space — selection clears

## Stacked bar filtering

1. Set Stack to SEX, **On Click** to **Filter**
2. Click on a stacked segment (e.g., Male within a RACE bar) — only matching rows remain
3. **Ctrl + click** another segment — that category is added to the filter
4. Click empty space — filter resets, all data returns

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

## Selected rows indicator

1. Set Split to RACE, **On Click** to **Select**
2. Open Context Panel > **Selection > Show Selected Rows** (enabled by default)
3. Select some rows in the grid — bars should show highlighted portions for selected rows
4. Disable **Show Selected Rows** — highlighted portions disappear
5. Re-enable — selected portions reappear

## Clipped bar indicators

1. Set Split to RACE, Value to AGE
2. Set custom **Value > Min** to a value that clips some bars (e.g., 40)
3. Verify wavy break and arrow indicators appear on clipped bars
4. Disable **Style > Show Clipped Bar Indicators** — indicators disappear
5. Re-enable — indicators reappear
6. Clear custom Min — bars return to normal

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

## Context menu regions

1. Right-click on the **value axis** — axis-specific menu appears (axis type, include nulls, axis font)
2. Right-click on the **chart area** — full context menu with Reset View, Orientation, Data, Order, Controls, Selection
3. Set Stack to SEX, right-click on the **legend** — legend-specific menu appears

## Row source

1. Open Context Panel > **Data > Row Source** — set to **All** — chart shows all rows regardless of filter
2. Open Filter Panel and filter some rows — chart still shows all rows
3. Set Row Source to **Filtered** — chart updates to show only filtered rows
4. Set Row Source to **Selected** — chart shows only selected rows

---
{
  "order": 100,
  "datasets": ["System:DemoFiles/demog.csv"]
}
