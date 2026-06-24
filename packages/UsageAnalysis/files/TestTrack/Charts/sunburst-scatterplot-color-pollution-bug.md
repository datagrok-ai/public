---
feature: charts
sub_features_covered:
  - charts.sunburst
  - charts.sunburst.inherit-from-grid
target_layer: playwright
coverage_type: edge
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

Bug-focused regression scenario for github-3412: opening the color editor
on a Scatterplot legend item removes color from a Sunburst that shares
the same color column. The pollution serializes into shared storage and
persists across viewer reopen. Fix landed in Charts 1.25 — this scenario
locks the cross-viewer color-state isolation invariant.

`pyramid_layer: bug-focused`. The canonical `sunburst.md` is single-viewer;
this scenario adds the cross-viewer pairing that is the bug-class
invariant. Color-editor popup driving requires DOM selectors deferred
per cycle charts-migrate-2026-05-07; programmatic substitute via
`grok.shell` color category APIs where available.

`related_bugs: [github-3412]` — bug-library reproduction class:
two viewers (Sunburst, Scatterplot) sharing a color column must not
pollute each other's color state when one's editor is opened.

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

- **github-3412 invariant carrier:** Sunburst's color binding survives
  Scatterplot's color-state mutation (Step 6 → Step 7). If the color
  editor DOM mutation is selector-pending, the programmatic substitute
  via `dataFrame.col.meta.colors` covers the bug-class invariant in
  effect (the bug is about shared storage corruption, which is
  observable via API).
- **No molecule data on SPGI.csv** — Stereo Category may not be
  directly present in SPGI.csv; the spec picks the first string column
  that fits the cross-viewer color pairing pattern.
- **Authority:** atlas-driven; closes the bug coverage gap for
  github-3412 surfaced in chain rev 2.

## Dataset metadata

```json
{
  "order": 36,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
```
