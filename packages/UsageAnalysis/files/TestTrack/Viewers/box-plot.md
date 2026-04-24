# Box plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Box plot

## Plot style: box vs violin

1. On the Context Panel > Value set Value to AGE
2. On the Context Panel > Data set Category 1 to RACE
3. On the Context Panel > Style and switch  Plot Style to `violin`
4. On the Context Panel > Style set:
  * **Bins** to 50
  * **Bins** to 500
  * **Interquartile Line Width** to 10
  * Switch Plot Style back to `box`

## Two-level categories

1. On the Context Panel > DataSet **Category 1** to SEX
2. On the Context Panel > DataSet **Category 2** to RACE
3. Toggle **Show Minor Categories** off
4. Toggle **Show Minor Categories** on
5. Enable **Show All Categories**
6. Disable **Show All Categories**
7. Set Category 2 to None

## Statistics display

1. Open Context Panel > **Statistics check Show Statistics** is enabled by default
2. Enable **Show Total Count**
3. Enable **Show Inliers Count**
4. Enable **Show Outliers Count**
5. Enable **Show Stdev**
6. Enable **Show Q1**
7. Enable **Show Q3**
8. Disable **Show Statistics** -- entire stats table disappears

## Markers

1. On the Context Panel > Value set Value to AGE
2. On the Context Panel > Data set Category 1 to SEX
3. On the Context Panel > Markers:
 - Change **Marker Type** to Square
 - Set **Marker Size** to 10
 - Set **Marker Opacity** to 80
 - Set **Markers Column** to RACE
 - Set **Marker Size Column** to WEIGHT
 - Set Markers Column to None
 - Set Marker Size Column to None

## Marker and Bin color coding

1. On the Context Panel > Color:
  - set Marker Color to AGE
  - **Color Axis Type** to `logarithmic`
  - Toggle **Invert Color Scheme** on
  - Set **Marker Color Column** to RACE -- categorical coloring applied
  - Toggle **Invert Color Scheme** off
  - Set Marker Color Column to None
  - Set Bin Color Column to WEIGHT
  - Set **Bin Color Aggr Type** to **min**
  - Set **Bin Color Aggr Type** to **max**
  - Set **Bin Color Aggr Type** to **med**
  - Set Bin Color Column to None

## Value axis configuration

1. On the Context Panel > **Value** set **Value** to AGE
2. On the Context Panel > **Value**:
   - Set **Axis Type** to **log**
   - Toggle **Invert Y Axis** on
   - Set **Min** to 20
   - Set **Max** to 60
   - Clear **Min** and **Max**
   - Set **Axis Type** back to **linear**
   - Toggle **Invert Y Axis** off

## Zoom by filter with filter panel

1. On the Context Panel > **Value** set **Value** to AGE
2. On the Context Panel > **Value** verify **Zoom Values By Filter** is enabled by default
3. Disable **Zoom Values By Filter**
4. Enable **Zoom Values By Filter**
5. On the Context Panel > **Data** set **Category 1** to RACE
6. Open Filter Panel
7. Narrow the AGE filter range
8. Verify viewer Y axis zoomed to match filter
9. Reset filter and close the Filter Panel

## Show empty categories

1. On the Context Panel > **Data** set **Category 1** to RACE
2. On the Context Panel > **Value** disable **Show Empty Categories**
3. Enable **Show Empty Categories**

## Box plot components

1. On the Context Panel toggle **Show Mean Cross** off
2. Toggle **Show Median Dash** off
3. Toggle **Show Upper Dash** off
4. Toggle **Show Lower Dash** off
5. Toggle **Show Inside Values** off
6. Toggle **Show Outside Values** off
7. Re-enable all six toggles

## Controls visibility

1. On the Context Panel > **Value** disable **Show Value Selector**
2. On the Context Panel > **Category** disable **Show Category Selector**
3. On the Context Panel > **Color** disable **Show Color Selector**
4. On the Context Panel > **Size** disable **Show Size Selector**
5. On the Context Panel > **Value** disable **Show Value Axis**
6. On the Context Panel > **Category** disable **Show Category Axis**
7. Re-enable all six controls

## Title and description

1. On the Context Panel > **Description** enable **Show Title**
2. Set **Title** to "Age by Race"
3. Set **Description** to "Box plot of patient ages"
4. Set **Description Visibility Mode** to Always
5. Change **Description Position** to Bottom
6. Set **Description Visibility Mode** to Never

## Date category mapping

1. On the Context Panel > **Data** set **Category 1** to STARTED
2. On the Context Panel > **Data** set **Category 1 Map** to **Month**
3. Set **Category 1 Map** to **Quarter**
4. Set **Category 1** to RACE

## Style customization

1. On the Context Panel > **Style**:
   - Set **Whisker Line Width** to 4
   - Set **Whisker Width Ratio** to 1.0
   - Set **Whisker Width Ratio** to 0.3
   - Toggle **Auto Layout** off
   - Toggle **Axis Use Column Format** off
   - Toggle **Auto Layout** on
   - Toggle **Axis Use Column Format** on

## Viewer filter formula

1. On the Context Panel > **Data** open **Filter**
2. Set filter formula to `${AGE} > 40`
3. Clear filter formula

## P-value (t-test)

1. On the Context Panel > **Data** set **Category 1** to SEX (2 categories)
2. On the Context Panel > **Value** set **Value** to AGE
3. On the Context Panel toggle **Show P Value** off
4. Toggle **Show P Value** on -- p-value should reappear
5. Press **T** key -- p-value should toggle off
6. Press **T** key -- p-value should toggle on
7. On the Context Panel > **Data** set **Category 1** to RACE (more than 2 categories) -- Alexander-Govern test should be used

## Legend

1. On the Context Panel > **Data** set **Category 1** to SEX
2. On the Context Panel > **Color** set **Marker Color Column** to RACE -- legend should show RACE categories
3. On the Context Panel > **Legend**:
   - Set **Legend Visibility** to **Never**
   - Set **Legend Visibility** to **Always**
   - Change **Legend Position** to Top Right
   - Change **Legend Position** to Bottom Left
4. On the Context Panel > **Color** set **Marker Color Column** to None

## Layout save and restore

1. On the Context Panel > **Value** set **Value** to WEIGHT
2. On the Context Panel > **Data** set **Category 1** to RACE
3. On the Context Panel > **Color** set **Marker Color Column** to SEX
4. On the Context Panel > **Statistics** enable **Show Total Count**
5. On the Context Panel > **Style** set **Plot Style** to violin
6. Save layout via JS API
7. Close all, reopen demog, add Box plot
8. Load saved layout via JS API -- verify Value is WEIGHT, Category 1 is RACE, Plot Style is violin, Marker Color Column is SEX, Show Total Count is enabled
9. Delete saved layout

## Visualization zoom (project and layout)

1. Set **Value** to AGE
2. Set **Category 1** to RACE
3. Zoom in
4. Apply **Reset View** -- zoom should reset
5. Zoom in
6. Save current project via JS API
7. Close all, reopen saved project
8. Verify that the zoom state is preserved in the Box plot viewer
9. Save layout via JS API
10. Close the viewer and load saved layout
11. Verify that the zoom state is NOT preserved in the layout
12. Delete saved layout and project

## Data properties and table switching (spgi-100)

1. Open spgi-100.csv, open demog
2. Add a box plot
3. On the Context Panel > **Data** set **Table** to spgi-100
4. On the Context Panel > **Value** set **Value** to Average Mass
5. On the Context Panel > **Data** set **Category 1** to Series
6. On the Context Panel > **Data** open **Filter** and set formula to `${Average Mass} > 225`
7. On the Context Panel > **Color** set **Bin Color Column** to TPSA
8. Close all

---
{
  "order": 8,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
