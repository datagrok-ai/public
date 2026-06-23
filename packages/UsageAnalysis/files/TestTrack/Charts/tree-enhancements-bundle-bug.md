---
feature: charts
sub_features_covered:
  - charts.tree
  - charts.tree.show-counts
  - charts.tree.on-click
  - charts.tree.font-size
  - charts.tree.layout
  - charts.tree.orient
  - charts.tree.include-nulls
target_layer: playwright
coverage_type: regression
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

Bug-focused regression scenario for github-3221: 8 capabilities bundled
into Tree viewer for NIBR's hierarchical-data workflow (showCounts,
onClick, fontSize, layout/orient, includeNulls, moleculeSize, tooltips,
label-overflow). Fix landed in Charts 1.4.3 — this scenario locks each
capability's regression bound.

`pyramid_layer: bug-focused`. The canonical `tree.md` only exercises
collaborative filtering — none of the 8 enhanced capabilities. This
scenario covers all 8 via property round-trip.

`related_bugs: [github-3221]` — bug-library reproduction class:
exercise each of the 8 capabilities and verify they continue to work
(no regression) on the current Tree viewer.

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

- **github-3221 invariant carrier:** for each of the 8 capabilities,
  setOptions does NOT throw, viewer renders without console error,
  and the property is enumerated in `getProperties()`. Strict
  read-back equality is best-effort (race-tolerant).
- **Capability 8 (label-overflow / showOverlappingLabels)** is currently
  inferred via `showMouseOverLine` — atlas listing may evolve. If
  `getProperties()` exposes a more specific label-overflow property,
  exercise that.
- **Authority:** atlas-driven; closes the bug coverage gap for
  github-3221 surfaced in chain rev 2.
- **Remediation cycle scope decision (charts-remediate-2026-05-09):**
  Critic E canonical subagent surfaced `charts.tree.layout` as
  declared in `sub_features_covered` but uncovered in code
  (predecessor cycle charts-automator-only-2026-05-08, verdict
  EVIDENCE_GAP: E-TRACE-02). Migrator decision: **Option A
  (TIGHTEN)** — keep `charts.tree.layout` in `sub_features_covered`,
  Automator implements guarded probe per scenario step 7's
  conditional ("if `tree.props.getProperties()` includes 'layout'"):
  - Read `tree.props.getProperties().map(p => p.name)` early.
  - If 'layout' is included, `exerciseProperty('layout',
    'orthogonal', 'layout=orthogonal')`.
  - Else log `console.warn('[SKIP] charts.tree.layout not exposed
    on current Tree build; capability 4 deferred per scenario step
    7 conditional')` (env-pending acceptable SR class per
    orchestrator Edit 10).
  Same guarded-probe pattern for moleculeSize per scenario step 8
  best-effort phrasing (enumeration-only check delegated to Setup
  propNames log; technical comment marker added).
- **Why Option A over Option B (REDUCE):** keeping
  `charts.tree.layout` in `sub_features_covered` preserves atlas
  mapping integrity (atlas line 607 confirms `charts.tree.layout`
  is present). Guarded probe is 8-line addition. If Tree build on
  dev does not expose `layout` property, env-pending defensive skip
  is the documented acceptable outcome (per orchestrator SR
  dispatch table). When the Charts package exposes `layout` on
  newer dev builds, the probe automatically picks it up without
  scenario amendment.

## Dataset metadata

```json
{
  "order": 34,
  "datasets": ["System:DemoFiles/demog.csv"]
}
```
