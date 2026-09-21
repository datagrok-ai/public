---
feature: piechart
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Automated coverage lives in the other scenarios in this folder; this file
  lists only checks needing a human eye.
---

# Pie chart tests (manual)

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
2. On the Filter Panel, tick only the **M** category in the SEX filter, then drag the AGE
   range down to roughly 40–60 — after each change the pie slices redraw to reflect only the
   filtered rows (slice proportions update immediately)

## Data section properties

1. Go to the Context Panel > Data section and exercise each property with a concrete change:
   set **Category** to SEX and back to RACE — the slices regroup by the chosen column; set
   **Row Source** to Selected and select a few grid rows — the chart rebuilds from the
   selection only; set **Filter** to `${AGE} > 40` — the slices reflect only matching
   rows; exercise the remaining Data properties the same way (e.g. the segment value
   column and its aggregation)

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv"]
}
