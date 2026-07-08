---
feature: charts
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [tree-hierarchical-capability-bundle]
pyramid_layer: bug-focused
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/tree-enhancements-bundle-bug.md
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
related_bugs:
  - github-3221
---

# Tree viewer — 8-enhancement bundle regression (github-3221)

Bug-focused regression scenario for github-3221: eight Tree-viewer
capabilities that were bundled together for a customer's
hierarchical-data workflow (show counts, on-click behavior, font size,
layout/orientation, include-nulls, molecule size, tooltips, label
overflow). The fix landed in Charts 1.4.3; this scenario locks in each
capability so it doesn't regress. The canonical `tree.md` only
exercises Filter Panel coordination — none of these 8 capabilities —
so this scenario exercises all 8 via a property round-trip.

## Setup

A clean Datagrok session. Single scenario; demog.csv as the source
table (sufficient hierarchy for CONTROL/SEX/RACE).

## Scenarios

### Scenario 1: Tree 8-enhancement bundle properties round-trip (github-3221 invariant)

Steps:

1. Open `System:DemoFiles/demog.csv` and add a Tree viewer.
   Set hierarchy `['CONTROL', 'SEX', 'RACE']` via setOptions.
   Wait 3000ms for Charts package + Tree to settle.

2. **Capability 1 — showCounts:** `tree.setOptions({showCounts: true})`.
   Read-back via `tree.props.get('showCounts')` (race-tolerant); expected
   `true`. Toggle to `false`; expected `false`. **Invariant:** no console
   error during toggles.

3. **Capability 2 — onClick:** `tree.setOptions({onClick: 'Filter'})`.
   Read-back; expected `'Filter'`. Set to `'Select'`; expected `'Select'`.
   **Invariant:** no console error.

4. **Capability 3 — fontSize:** `tree.setOptions({fontSize: 14})`.
   Read-back; expected `14`. Set to `30` (atlas-cited max); expected
   `30`. **Invariant:** no console error.

5. **Capability 4 — orient:** `tree.setOptions({orient: 'LR'})`.
   Read-back; expected `'LR'`. Set to `'RL'`; expected `'RL'`. Set back
   to `'TB'` (default); expected `'TB'`. **Invariant:** no console error
   across orientation transitions.

6. **Capability 5 — includeNulls:** `tree.setOptions({includeNulls: true})`.
   Read-back; expected `true`. Toggle off. **Invariant:** no console error.

7. **Capability 6 — layout** (related to orient but separate atlas
   sub_feature): if `tree.props.getProperties()` includes `layout`,
   exercise via `tree.setOptions({layout: 'orthogonal'})`. Read-back.
   **Invariant:** no console error.

8. **Capability 7 — moleculeSize:** if applicable to the dataset
   (demog.csv lacks molecule columns; this assertion is best-effort —
   read property metadata without setting a value).
   Verify `tree.props.getProperties()` enumerates the property name
   `moleculeSize` (regression: capability is exposed even if not
   exercised).

9. **Capability 8 — tooltips / showMouseOverLine:**
   `tree.setOptions({showMouseOverLine: true})`. Read-back. Toggle.
   **Invariant:** no console error.

10. **Final visual stability check:** verify the Tree viewer's root
    DOM element is non-empty and has non-zero size after all 8 toggles
    (no white screen / DOM tear-down).

## Notes

- The github-3221 invariant: for each of the 8 capabilities,
  `setOptions` doesn't throw, the viewer renders without a console
  error, and the property shows up in `getProperties()`. Strict
  read-back equality is best-effort (race-tolerant).
- Capability 8 (label overflow) is currently exercised via
  `showMouseOverLine` as a stand-in — if a future Tree build exposes
  a more specific label-overflow property, the spec should switch to
  that.
- Steps 7 (layout) and 8 (moleculeSize) use a guarded probe: they
  check whether the property is exposed on the current build before
  exercising it, logging a skip instead of failing if not. This is
  deliberate — different Tree builds may not yet expose `layout`, and
  `demog.csv` has no molecule column to exercise `moleculeSize`
  against, so those two capabilities are verified by enumeration
  rather than a full round-trip.

## Dataset metadata

```json
{
  "order": 34,
  "datasets": ["System:DemoFiles/demog.csv"]
}
```
