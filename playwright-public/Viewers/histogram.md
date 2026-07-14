# Histogram tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Histogram

## Bins configuration

1. Set Value to AGE
2. Open Context Panel > **Value > Bins** -- default is 20
3. Set Bins to 5
4. Set Bins to 100
5. Set Bins to 1
6. Reset Bins to 20
7. Set **Bin Width Ratio** to 1.0
8. Set **Bin Width Ratio** to 0.3
9. Set **Bin Width Ratio** to 0.8

## Split column

1. Open Context Panel > **Category > Split Column** > set to SEX
2. Enable **Normalize Values**
3. Disable **Normalize Values**
4. Toggle **Show Markers** off
5. Set **Spline Tension** to 5
6. Set Split Column to RACE
7. Set Split Column to None
8. Set Split Column to SEX
9. Enable **Split Stack**
10. Enable **Value > Show Values**
11. Disable **Split Stack**
12. Enable **Misc > Show Distribution Lines**
13. Disable **Show Distribution Lines**

## Color coding

1. Set Split Column to None
2. Open Context Panel > **Color > Color Column** > set to WEIGHT
3. Set **Color Aggr Type** to min
4. Set **Color Aggr Type** to max
5. Toggle **Invert Color Scheme** on
6. Toggle **Invert Color Scheme** off
7. Set Color Column to None

## Value range

1. Set Value to AGE
2. Set **Value > Min** to 30
3. Set **Value > Max** to 60
4. Clear Min and Max
5. Enable **Show Range Inputs**
6. Type 30 in the min input and 60 in the max input
7. Disable **Show Range Inputs**

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

## Bin selection

1. Set Value to AGE
2. Click a bin -- rows in that bin should be selected, status bar shows selected count
3. Ctrl+click another bin -- selection of first bin should invert
4. Click empty area -- selection persists
5. Shift+drag across several bins -- all bins in the drag range should be selected
6. Click a single bin again -- only that bin's rows are selected
7. Set Split Column to SEX
8. Enable **Split Stack**
9. Click on a specific color segment in a stacked bin -- only rows matching that bin AND category should be selected
10. Disable **Split Stack**
11. Click on a spline line -- rows for that category should be selected

## Filtering

1. Apply a filter on SEX column to "M" only via JS API
2. Toggle **Data > Show Filtered Out Rows** off -- only filtered bars remain
3. Toggle **Show Filtered Out Rows** on -- semi-transparent bars return
4. Change **Style > Filtered Out Color** -- filtered-out bar color should update
5. Remove the filter
6. Set Value to AGE
7. Go to the Context Panel > Data, enable Filtering Enabled
8. Drag the horizontal range slider to narrow the range — verify the dataset filters accordingly (row count in status bar decreases)
9. Toggle **Normalize To Filter** off -- bar heights should reflect global counts
10. Toggle **Normalize To Filter** on -- heights normalize to current filter
11. Toggle **Zoom To Range** off -- X axis shows full column range even when filtered
12. Toggle **Zoom To Range** on -- X axis zooms to visible range
13. Enable **Bin To Range** -- bins recalculate based on filtered range only
14. Disable Filtering Enabled — verify the dataset is no longer filtered (row count restores)
15. Open the Filter Panel and check Asian and Others categories in the RACE filter — verify the histogram updates
16. Close the Filter Panel

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
