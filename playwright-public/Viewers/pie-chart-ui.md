# Pie chart tests (manual)

Ручной чеклист. Не входит в автоматизацию PW.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a Pie chart

## Tooltips

1. Hover over a slice — tooltip appears showing category name, row count, and column statistics
2. Set Segment Angle Column to AGE — hover over a slice, tooltip should include AGE aggregation
3. Set Segment Length Column to WEIGHT — tooltip should also show WEIGHT aggregation
4. Move mouse away from pie — tooltip disappears

## Filtering (exploratory)

1. Go to Tables and click the Filter icon
2. Change arbitrary settings on the Filter Panel and check the pie chart interaction

## Check all data properties

1. Go to the Context Panel > Data section and check all the properties

## Range slider

1. Add a title and a description
2. Check the range slider functionality

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv"]
}
