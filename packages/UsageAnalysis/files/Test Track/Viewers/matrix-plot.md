# Matrix plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Matrix plot

## Default State

1. Verify the Matrix plot viewer is present in the view.
2. Verify default **X Column Names** contains only numerical columns (up to 10).
3. Verify default **Y Column Names** contains only numerical columns (up to 10).

## Column Configuration

1. Set **X Column Names** to AGE and HEIGHT via JS API.
2. Verify the matrix shows only AGE and HEIGHT on the X axis.
3. Set **Y Column Names** to AGE, HEIGHT, WEIGHT via JS API.
4. Verify the matrix shows a 2×3 grid of cells.
5. Reset **X Column Names** and **Y Column Names** to defaults via JS API.

## Cell Plot Type

1. Open the property panel (click the gear icon on the Matrix plot title bar).
2. Set **Cell Plot Type** to "Scatter plot".
3. Set **Cell Plot Type** back to "Density plot".

## Show Axes

1. Open the property panel.
2. Set **Show X Axes** to true.
3. Set **Show Y Axes** to true.
4. Set **Show X Axes** to false.
5. Set **Show Y Axes** to false.

## Auto Layout

1. Open the property panel.
2. Verify **Auto Layout** is enabled by default.
3. Set **Show X Axes** and **Show Y Axes** to true.
4. Set **Auto Layout** to false.
5. Verify axes remain visible.
6. Set **Auto Layout** back to true.

## Scrolling

1. Set **X Column Names** to AGE, HEIGHT, WEIGHT, STARTED via JS API.
2. Set **Y Column Names** to AGE, HEIGHT, WEIGHT, STARTED via JS API.
3. Verify horizontal and vertical scrollbars are visible in the viewer.
4. Drag the horizontal scrollbar handle to the right.
5. Verify column labels on the X axis change to reflect the new viewport position.
6. Drag the vertical scrollbar handle downward.
7. Verify column labels on the Y axis change to reflect the new viewport position.

## Row Source

1. Open the property panel.
2. Change **Row Source** to "Selected".
3. In the grid, select 50 rows via JS API.
4. Verify the matrix plot updates to reflect only the selected rows.
5. Change **Row Source** back to "Filtered".

## Filter Integration

1. Open the filter panel.
2. In the SEX filter, select only "M".
3. Verify the filter count decreases.
4. Remove the SEX filter.
5. Verify the filter count returns to the total row count.

## Data Filter Property

1. Open the property panel.
2. In the **Filter** field (Data category), enter `${AGE} > 30`.
3. Verify the row count shown in the plot decreases.
4. Clear the **Filter** field.
5. Verify the plot returns to showing all rows.

## Back Color

1. Set **Back Color** to red via JS API.
2. Verify the background color of the Matrix plot changes.
3. Reset **Back Color** to white via JS API.

## Title and Description

1. Open the property panel.
2. Set **Show Title** to true.
3. Set **Title** to "My Matrix".
4. Verify the title "My Matrix" appears on the viewer.
5. Set **Description** to "Test description".
6. Set **Show Title** to false.

## Inner Viewer Look

1. Open the property panel. Set **Cell Plot Type** to "Scatter plot".
2. Set inner viewer marker size to 10 via JS API.
3. Verify all off-diagonal scatter plot cells update with the new marker size.
4. Set **Cell Plot Type** back to "Density plot".
5. Set inner viewer bin count to 20 via JS API.
6. Verify all off-diagonal density plot cells update accordingly.

## Context Menu

1. Right-click on the Matrix plot viewer.
2. Select **Clone** from the context menu.
3. Verify a copy of the Matrix plot is added to the view.
4. Close the cloned viewer.

---
{
  "order": 16,
  "datasets": ["System:DemoFiles/demog.csv"]
}
