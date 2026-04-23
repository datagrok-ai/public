# Heat map — manual test checklist

Ручной чеклист. Не входит в автоматизацию PW.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Heat map

## Alt+drag area zoom

1. Alt+drag over a rectangular area on the heat map — view zooms into that area
2. Verify horizontal and vertical range sliders update to reflect the zoomed region
3. Alt+drag over a smaller area within the zoomed view — zooms in further
4. Double-click range sliders to reset the view to full extent

## Right-click content panning

1. Zoom in using Alt+drag, then right-click and drag on the content area — content pans
2. Release right-click — panning stops, content stays at new position

## Custom sorting on categorical column

> Note: requires SPGI dataset.
> Setup: Close all, open SPGI.csv, go to Tables > SPGI, add Heat map.

1. Right-click the column header Primary Series Name
2. Select Sort > Custom, move the empty value to the top of the list, click OK
3. Right-click the same column, choose Sort > Ascending — empty value appears first
4. Right-click the same column, choose Sort > Descending — sort order is reversed

---
{
  "order": 14,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
