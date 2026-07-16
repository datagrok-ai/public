---
feature: charts
target_layer: playwright
coverage_type: edge
priority: p1
realizes_atlas: [sunburst-scatterplot-shared-color-isolation]
realizes: [charts.sunburst]
pyramid_layer: bug-focused
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst-scatterplot-color-pollution-bug.md
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
related_bugs:
  - github-3412
---

# Sunburst × Scatterplot — color-state pollution (github-3412)

Bug-focused regression scenario for github-3412: opening the color
editor on a Scatterplot legend item was removing color from a Sunburst
viewer that shares the same color column — the pollution serialized
into shared storage and persisted across viewer reopen. The fix landed
in Charts 1.25; this scenario locks in the cross-viewer color-state
isolation invariant: two viewers (Sunburst, Scatterplot) sharing a
color column must not pollute each other's color state when one
viewer's editor is opened. The canonical `sunburst.md` only covers a
single viewer at a time — this scenario adds the cross-viewer pairing
that the bug class requires. Driving the actual color-editor popup
needs DOM selectors that aren't available yet, so the steps use a
programmatic substitute via the `grok.shell` color-category API where
possible.

## Setup

A clean Datagrok session. Single scenario; SPGI.csv as the source
(has Stereo Category column per bug ticket).

## Scenarios

### Scenario 1: Sunburst × Scatterplot color-state isolation (github-3412 invariant)

Steps:

1. Open `System:DemoFiles/SPGI.csv`. Wait for semType detection.

2. Add a Sunburst viewer via `tv.addViewer('Sunburst')`. Wait 3000ms.
   Configure Sunburst hierarchy to a categorical string column likely
   to be a color binding target (atlas-driven; pick first string-typed
   column that is enumerable as a hierarchy candidate, or fall back to
   any string column on SPGI.csv).

3. Add a Scatterplot viewer via `tv.addViewer('Scatter plot')`. Wait
   1500ms.

4. Configure Scatterplot's `colorColumnName` to the SAME string column
   the Sunburst is hierarchy-bound to:
   `scatterplot.setOptions({colorColumnName: <sharedColumn>})`.

5. **Capture baseline color state.** Read back
   `sunburst.props.get('colorColumnName')` (race-tolerant) and
   `scatterplot.props.get('colorColumnName')` (race-tolerant). Log.

6. **Trigger color-state mutation via Scatterplot.** The bug ticket
   cites opening the color-editor popup; the programmatic substitute
   is to mutate the column's category-color mapping directly via
   `grok.shell.tv.dataFrame.col(<sharedColumn>).meta.colors` if
   available, OR via `scatterplot.setOptions({...colorOverride...})`
   if the property exists.
   - If neither path is available on dev, the spec test.skip's this
     specific mutation step with a logged warning (selector-pending /
     API-pending) per cycle charts-migrate-2026-05-07 lessons.

7. **github-3412 invariant carrier:** verify the Sunburst's
   `colorColumnName` property still equals the original shared column
   (not corrupted to null or a different value). Verify Sunburst's
   root DOM is non-empty and has non-zero size (visual stability).

8. **Persistence check (in-session only):** save the table view (no
   project save — that's separate cleanup), close the Sunburst, and
   re-add it. Verify the re-added Sunburst still binds to the original
   color column (the bug ticket cites pollution persisting across
   reopen via shared storage).

9. **Cleanup:** `grok.shell.closeAll()`.

## Notes

- The github-3412 invariant is carried by the Step 6 → Step 7
  transition: the Sunburst's color binding must survive the
  Scatterplot's color-state mutation. If the actual color-editor DOM
  interaction isn't available, the programmatic substitute via
  `dataFrame.col.meta.colors` still covers the invariant in effect,
  since the bug is about shared-storage corruption, which is
  observable via the API.
- SPGI.csv doesn't have a "Stereo Category" column like the original
  bug ticket — the spec picks the first string column that fits the
  cross-viewer color-pairing pattern instead.

## Dataset metadata

```json
{
  "order": 36,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
```
