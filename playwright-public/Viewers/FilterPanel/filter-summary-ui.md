---
feature: filters
realizes_atlas: [filters.int.active-counter-counts-filtering-only]
realizes: []
priority: p1
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  The filter-summary indicator on the Filter Panel header is a visual
  affordance whose current form (a numeric badge; formerly a question-mark
  icon) needs visual confirmation of the listed content.
related_bugs: []
---

# Filter Panel — filter summary indicator

Filtering and Reset are covered automatically; this file covers only
the summary indicator content.

## Steps

1. Close all, open spgi-100
2. Open the Filter Panel and apply two filters: Stereo Category (categorical)
   and Average Mass (numerical range)
3. Add a scatterplot, barchart, pie chart, trellis plot, PC plot, Scaffold tree and use them to filter the dataset
4. Check the filter-summary indicator on the Filter Panel header (numeric
   badge) — it reflects all active filters; hover it — the details list
   all filters with their categories/range
5. Click the Reset filter icon — the indicator empties and all filtering is
   reset (full row count restored)
---
{
  "order": 17,
  "datasets": ["System:AppData/Chem/tests/spgi-100.csv"]
}