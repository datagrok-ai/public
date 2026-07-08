# 3D Scatter Plot tests — Playwright

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add 3D Scatter Plot

## Axis column assignment

1. Open properties
2. Set X to AGE, Y to HEIGHT, Z to WEIGHT
3. Set X to WEIGHT, Y to AGE, Z to HEIGHT
4. Set X back to AGE, Y to HEIGHT, Z to WEIGHT

## Axis types

1. Open properties, set X Axis Type to logarithmic
2. Set Y Axis Type to logarithmic
3. Set Z Axis Type to logarithmic
4. Set X Axis Type back to linear
5. Set Y Axis Type back to linear
6. Set Z Axis Type back to linear

## Color coding — categorical

1. Open properties, set Color to SEX
2. Set Color to RACE
3. Clear Color

## Color coding — numerical

1. Open properties, set Color to AGE
2. Clear Color

## Size coding

1. Open properties, set Size to WEIGHT
2. Set Size to AGE
3. Clear Size

## Labels

1. Open properties, set Label to SEX
2. Clear Label

## Marker type

1. Open properties, set Marker Type to sphere
2. Set Marker Type to box
3. Set Marker Type to cylinder
4. Set Marker Type to tetrahedron
5. Set Marker Type to dodecahedron
6. Set Marker Type back to octahedron

## Marker opacity and rotation

1. Open properties, set Marker Opacity to 20
2. Set Marker Opacity to 100
3. Set Marker Opacity back to 69
4. Enable Marker Random Rotation
5. Disable Marker Random Rotation

## Filtered out points

1. Open the filter panel
2. Filter AGE to range 20–40
3. Open properties, enable Show Filtered Out Points
4. Disable Show Filtered Out Points
5. Clear the AGE filter

## Axes visibility and grid lines

1. Open properties, disable Show Axes
2. Enable Show Axes
3. Disable Show Vertical Grid Lines
4. Disable Show Horizontal Grid Lines
5. Enable Show Vertical Grid Lines
6. Enable Show Horizontal Grid Lines

## Background and colors

Note: color properties set via JS API (UI color picker is not automatable).
1. Set Back Color to black
2. Set Axis Line Color to white
3. Restore Back Color to white and Axis Line Color to default

## Dynamic camera movement

1. Open properties, enable Dynamic Camera Movement
2. Disable Dynamic Camera Movement

## Zoom and navigation

1. Scroll mouse wheel up over the plot five times
2. Scroll mouse wheel down five times
3. Right-click on the plot, click Reset View

## Mouse-over row group highlight

1. Add a Bar Chart to the view
2. Open 3D Scatter Plot properties, confirm Show Mouse Over Row Group is checked
3. Move mouse over a bar in the Bar Chart
4. Open 3D Scatter Plot properties, uncheck Show Mouse Over Row Group
5. Move mouse over a bar in the Bar Chart
6. Recheck Show Mouse Over Row Group
7. Close the Bar Chart

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
