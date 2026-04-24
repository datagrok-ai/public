# Legend color consistency — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI | 10s | PASS | PASSED | 3624 rows |
| 2 | Add 6 viewers (Hist, Line, Bar, Pie, Trellis, Box) — omit Scatter | 5s | PASS | PASSED | tv.viewers = [Grid, Hist, Line, Bar, Pie, Trellis, Box] |
| 3 | Set Stereo Category legend on each viewer | 4s | PASS | PASSED | Hist/Line/Bar→split, Pie→category, Trellis→X, Box→categoryColumnNames |
| 4 | Enable categorical color coding + change 2 colors (R_ONE→red, S_UNKN→green) | 3s | PARTIAL | PASSED | R_ONE correctly becomes red (0xffff0000); S_UNKN does NOT become green — `col.meta.colors.setCategorical({key:value})` iterates positionally so colors leak across categories |
| 5 | Pick new color in Bar chart legend color picker | 1s | SKIP | n/a | Bar chart legend DOM is absent even with Always visibility; same Box plot; color picker not reachable on these viewers |
| 6 | Save layout → re-apply — palette persists | 5s | PASS | PASSED | `.categorical-colors` tag and `meta.colors.getColor(0)` both survive round-trip |
| 7 | Save project, reopen | 2s | FAIL | PASSED (asserted against error) | `project_relations_entity_id_fkey` FK constraint — known limitation |
| 8 | Cleanup | 1s | PASS | n/a | Layout deleted, closeAll |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 55s |
| grok-browser execution (scenario steps) | 35s |
| Execute via grok-browser (total) | 2m 30s |
| Spec file generation | 35s |
| Spec script execution | 26s |
| **Total scenario run (with model)** | 3m 31s |

## Summary

Color consistency through layout round-trip works — the `.categorical-colors` tag survives save/reload and `R_ONE` stays red. Two gaps: `col.meta.colors.setCategorical({name: color})` does not honor the category-name keys (positional iteration produces wrong assignments for keys beyond the first), and Bar chart and Box plot don't surface a legend DOM even with `Always` visibility so color-picker interactions can't be verified. Project save fails with the known FK constraint. **Total scenario run (with model): 3m 31s**.

## Retrospective

### What worked well
- `.categorical-colors` tag applied directly on the column propagates
- Layout round-trip preserves the palette and `meta.colors.getColor()` readings

### What did not work
- `col.meta.colors.setCategorical(nameKeyedMap)` yields wrong mapping for categories after the first key (R_ONE correct, S_UNKN not green)
- Bar chart / Box plot have no DOM legend even with Always — color-picker step unreachable
- Project save: FK constraint (same as other scenarios)

### Suggestions for the platform
- Fix `ColumnColorHelper.setCategorical` to honor object keys by category name, not by positional iteration
- Make Bar chart and Box plot render a DOM legend when Always is set (currently the legend is absent)
- Auto-save unsaved child dataframes on `projects.save` or emit a typed error

### Suggestions for the scenario
- Document the preferred JS-API call for setting category colors (`col.meta.colors.setCategorical(array)` vs tag) so tests are reproducible
- Note that Bar chart legend color picker isn't reachable — suggest using a viewer that exposes a DOM legend (Histogram, Pie chart) for step 5
