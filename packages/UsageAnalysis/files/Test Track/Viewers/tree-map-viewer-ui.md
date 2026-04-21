Ручной чеклист. Не входит в автоматизацию PW.

# Tree Map — manual test checklist

## Hover tooltip

1. Open demog, add Tree Map
2. Set split column to RACE
3. Hover over a rectangle — a tooltip appears showing the leaf name and a sample of rows in that group
4. Move to a different rectangle — tooltip updates to the new group
5. Move off the canvas — tooltip disappears

## Default Color

1. Open demog, add Tree Map
2. Clear the split column so no category color applies
3. Open Settings → General → Default Color
4. Change Default Color to a custom value (e.g. light blue)
5. Verify all leaf rectangles render in the selected color

## Visual canvas rendering — color scale

1. Open demog, add Tree Map
2. Set split column to RACE
3. Set color column to AGE with aggregation avg
4. Verify rectangles are colored on a gray→red scale: groups with higher avg AGE appear more red
5. Change aggregation to max — verify the color distribution updates accordingly

## Visual canvas rendering — nested borders

1. Open demog, add Tree Map
2. Set first split column to RACE, add SEX as a second level
3. Verify inner borders separate the SEX sub-groups within each RACE rectangle
4. Verify outer borders delineate the top-level RACE groups
5. Verify leaf captions are centered inside each rectangle

## Tooltip content

1. Open demog, add Tree Map
2. Set split column to RACE, color column to AGE
3. Hover over a rectangle — tooltip header shows the leaf name
4. Verify tooltip body lists individual rows belonging to that group
