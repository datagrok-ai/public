# Density plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Density plot

## Zoom and Reset View

1. Zoom in with mouse wheel
2. Right-click the viewer and in the context menu click "Reset View"
3. Verify: full data range is visible again

## Axis Column Assignment

1. On the viewer, click the X column selector and set X column to WEIGHT
2. On the viewer, click the Y column selector and set Y column to HEIGHT
3. On the viewer, click the X column selector and set X column to AGE
4. On the viewer, click the Y column selector and set Y column to WEIGHT

## Bins

1. Go to the Context Panel > Misc and set:

 * Bin Shape to `rectangle`
 * Bin Shape to `hexagon`
 * Bins to 5
 * Bins to 100
 * Bins to 200
 * Bins to 50

## Show/Hide Color Scale

1. Go to the Context Panel > Style and set: Set Show Color Scale to false
2. Set Show Color Scale to true

## Axis Visibility

1. Go to the Context Panel > X and set Show X Axis to false
2. Go to the Context Panel > Y and set Show Y Axis to false
3. Go to the Context Panel > X and set Show X Axis to true
4. Go to the Context Panel > Y and set Show Y Axis to true

## Axis Inversion and Logarithmic Axis

1. Go to the Context Panel > X and set Invert X Axis to true
2. Set X Axis Type to `logarithmic`
3. Go to the Context Panel > Y and set Invert Y Axis to true
4. Set Y Axis Type to `logarithmic`
5. Go to the Context Panel > X and set Invert X Axis to false
6. Set X Axis Type to `linear`
7. Go to the Context Panel > Y and set Invert Y Axis to false
8. Set Y Axis Type to `linear`


## Show/Hide Selectors and Bin Slider

1. Go to the Context Panel > X and set Show X Selector to false
2. Go to the Context Panel > Y and set Show Y Selector to false
3. Go to the Context Panel > Misc and set Show Bin Selector to false
4. Go to the Context Panel > X and set Show X Selector to true
5. Go to the Context Panel > Y and set Show Y Selector to true
6.  Go to the Context Panel > Misc and set Show Bin Selector to true

## Min/Max Axis Bounds

1. Set X Min to 20
2. Set X Max to 50
3. Set Y Min to 100
4. Set Y Max to 200
5. Clear X Min and X Max
6. Clear Y Min and Y Max

## Color Transform Type

1. Set Color Transform Type to logarithmic
2. Set Color Transform Type to linear

## Title and Description

1. Enable Show Title
2. Set Title to "Density Distribution"
3. Set Description to "AGE vs HEIGHT density"
4. Set Description Visibility Mode to Always
5. Set Description Position to Bottom

## Selection

1. Click any bin in the viewer - Verify: some rows are selected 
2. Esc - Verify: selection is cleared

## Layout Persistence

1. Set X column to WEIGHT
2. Set Y column to HEIGHT
3. Set Bins to 25
4. Set Bin Shape to rectangle
5. Set Invert Color Scheme to true
6. Save the layout via JS API
7. Close Density plot viewer
8. Apply the saved layout
9. Verify density plot restores with WEIGHT, HEIGHT, 25 bins, rectangle shape, and inverted color scheme

## Row Source Filtering

1. Set Filter to `${AGE} > 30`
2. Clear Filter
3. Open spgi-100
4. Go back to the demog table and click the density plot
5. On the Context Panel > Data set Table to spgi-100 - check, no errors
6. Set table back to the demog
7. Close All
---
{
  "order": 11,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
