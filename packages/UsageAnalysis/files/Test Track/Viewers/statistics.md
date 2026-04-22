# Statistics viewer tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Statistics

## Add viewer

1. Click the Statistics icon in the Toolbox Viewers section — Statistics viewer opens
2. Close the Statistics viewer
3. Click the Statistics icon in the Toolbox Viewers section again — Statistics viewer opens

## Default statistics display

1. Verify the viewer shows stat columns: values, nulls, unique, min, max, avg, med, stdev
2. Verify numerical columns AGE, HEIGHT, WEIGHT have non-empty values in all stat columns
3. Verify the name column lists the demog dataset columns

## Statistics for categorical columns

1. Verify numerical stat columns (avg, min, max, stdev) are empty for categorical column rows (SEX, RACE, DIS_POP)
2. Verify count-type stat columns (values, nulls, unique) are populated for categorical column rows

## Add and remove statistics columns

1. Right-click the Statistics viewer
2. Hover over Statistics in the context menu — submenu opens listing available stat types
3. Click sum in the submenu — a sum column appears in the viewer
4. Right-click the viewer and hover over Statistics in the context menu
5. Click sum again — the sum column disappears

## Histogram columns

1. Right-click the Statistics viewer
2. Hover over Histograms in the context menu — submenu lists categorical columns with fewer than 10 categories (SEX, RACE, DIS_POP)
3. Click SEX in the submenu — a histogram sparkline column appears in the viewer
4. Right-click the viewer and hover over Histograms in the context menu
5. Click SEX again — the histogram column disappears

## Row source: filtered rows

1. Open the filter panel
2. Add an AGE range filter and set the range to 20–40
3. Verify the values count in the Statistics viewer decreases compared to the full dataset
4. Remove the AGE filter — statistics revert to full dataset values

## Row source: selected rows

1. Select several rows in the grid
2. Open the Properties panel for the Statistics viewer
3. Set Rows to Selected — statistics update to reflect only the selected rows
4. Verify the values count matches the number of selected rows
5. Set Rows back to All — statistics revert to full dataset

## Open as table

1. Right-click the Statistics viewer
2. Click Open as table — a new table tab opens with stat types as column headers and data columns as rows

## Full-screen shortcut

1. Click the expand icon on the Statistics viewer panel — viewer expands to full screen
2. Press Alt+F — viewer returns to normal size

---
{
  "order": 100,
  "datasets": ["System:DemoFiles/demog.csv"]
}
