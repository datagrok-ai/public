---
feature: scatterplot
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-17227
    status: fixed
  - GROK-19083
realized_as:
  - scatterplot-legend-spec.ts
  - missing-spec.ts
---

# Scatter plot legend

Legend entries follow the Color and Marker columns.

## Steps

1. Open demog
2. Add a scatter plot
3. Set Color to RACE

Expected: one legend entry per category.
