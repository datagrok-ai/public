# Pie chart tests (manual)

Automated coverage lives in the other scenarios in this folder; this file
lists only checks needing a human eye.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a Pie chart

## Tooltips

1. Hover over a slice — the slice is highlighted with the hover color, and the
   tooltip identifies that slice's category with its values

## Mouse-over cross-highlight

1. Set Category to RACE and turn on **Show Mouse Over Row Group**; add a Bar chart on the same table and hover one of its bars — the pie arc of the matching category is highlighted
2. Hover a pie slice — the rows contributing to that slice are highlighted in the grid and other viewers

## Filtering (exploratory)

1. Go to Tables and click the Filter icon
2. Change arbitrary settings on the Filter Panel and check the pie chart interaction

## Check all data properties

1. Go to the Context Panel > Data section and check all the properties

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv"]
}
