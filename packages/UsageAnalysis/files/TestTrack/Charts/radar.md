---
feature: charts
sub_features_covered:
  - charts.radar
  - charts.radar.show-current-row
  - charts.radar.color-column
  - charts.radar.color-palette
  - charts.echart-base.table
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-viewer-radar
  - viewer-property-panel-gear
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/radar.md
migration_date: 2026-05-07
migration_report: radar-migration-report.md
related_bugs:
  - GROK-18085
---

# Radar viewer (Charts package)

Smoke scenario for the Charts package Radar viewer. Verifies the
"Add Viewer → Radar" entry point on two distinct datasets and the
property panel surface (Gear icon → Context Panel) for the four
representative properties the manual scenario calls out: switching
the bound table, selection check-boxes, the chosen Values count, and
style/color changes.

`pyramid_layer: ui-smoke` per `scenario-chains/charts.yaml` rev 1
(Rule 1 — single viewer create + configure smoke; shortest qualifying
scenario by step count). This scenario owns the canonical
Add-Viewer + Property-Panel-Gear smoke surface for the Charts
section; sister scenarios `sunburst.md` and `tree.md` own their
specialty UI surfaces directly.

`related_bugs: [GROK-18085]` — Radar table-rebind on project
save/reopen. The original scenario does NOT exercise project
save/reopen, so the cross-cutting bug invariant is NOT verified by
this scenario; the citation surfaces awareness only. Coverage gap
flagged in chain rev 1 `bug_match_attempts_skipped` (skip_category:
reproduction_unparseable) and surfaced again in this migration's
report for Critic F downstream review.

## Setup

A clean Datagrok session is the only shared setup. Each step opens
its own dataset; no fixture chaining required (`depends_on: []` per
chain rev 1).

## Scenarios

### Add Radar viewer to two datasets and exercise the property panel

1. Open `System:DemoFiles/geo/earthquakes.csv` (e.g. via
   **Browse** > **Files**, or **File** > **Open**).
2. On the Menu Ribbon, click **Add viewer** and select **Radar**.
3. **Verify:** the Radar viewer opens for `earthquakes.csv` with
   no errors.
4. Open `System:DemoFiles/demog.csv`.
5. On the Menu Ribbon, click **Add viewer** and select **Radar**.
6. **Verify:** the Radar viewer opens for `demog.csv` with no errors.
7. On the Radar viewer (on `demog.csv`), click the **Gear** icon.
8. **Verify:** the **Context Panel** opens with the Radar viewer's
   properties.
9. **Switching tables:** in the Context Panel, change the bound
   `table` property between `earthquakes` and `demog`.
   **Verify:** the Radar viewer rebinds and re-renders against the
   new table without errors.
10. **Check-boxes in selection:** select two or more rows in the
    bound table (Ctrl+Click in the grid, or
    `df.selection.set(...)` equivalent UI gesture).
    **Verify:** the selected-row lines are reflected on the Radar
    viewer (per `charts.radar.show-current-row` family — current /
    mouseover row lines).
11. **Increasing and decreasing the amount of chosen Values:** in
    the Context Panel, change the columns chosen as Values (the
    color column / values bag).
    **Verify:** the Radar viewer redraws with more / fewer axes
    matching the chosen Values.
12. **Style (color) changes:** in the Context Panel, change a
    Style-category color (e.g. `currentRowColor`,
    `mouseOverRowColor`, `lineColor`, or `backgroundMinColor` /
    `backgroundMaxColor`).
    **Verify:** the Radar viewer redraws with the new color.
13. **Check all properties (broad sweep):** scroll through the
    Context Panel and confirm every Radar property is visible
    and editable. Spot-check that toggling representative
    properties does not throw console errors. (See Notes — this
    step preserves the original scenario's broad "Check all
    properties" instruction without enumerating each radar.*
    property; the four bullets above are the explicit Main things
    the original calls out as MUST be reflected on the viewer.)

## Notes

- **target_layer: playwright** — chosen because a sibling
  `radar-spec.ts` already exists at the playwright layer (per
  `existing-test-index.yaml` line 32099). Chain rev 1 deferred this
  decision to Step 2; selecting `playwright` aligns with the
  established sibling spec (rather than the chain's tentative
  `uitests-package` reading).
- **"Check all properties" (Step 13)** is preserved deliberately
  broad — the original scenario does not enumerate which of the 9
  atlas-listed `charts.radar.*` properties (title, min/max
  percentile, show-current-row, show-tooltip, color-column,
  color-palette, show-min-max, legend-visibility) are required and
  which are optional. Migrator preserves the original ambiguity
  rather than narrowing or widening (chain rev 1 unresolved
  ambiguity #1; flagged for QA review).
- **GROK-18085 (related_bug, not exercised here):** the bug
  reproduction requires "save the project, then reopen it" with a
  Radar viewer's bound table changed mid-session. The original
  radar.md has no project save / reopen step (3 numbered steps
  cover Add Viewer × 2 + property panel only). This scenario
  surfaces the cross-cutting bug citation but does NOT verify the
  invariant; chain rev 1 records this as a coverage gap
  (bug_match_attempts_skipped: reproduction_unparseable) for
  Critic F surfacing.
- **Helpers (existing in registry, available for downstream
  Automator):** `softStep`, `loginToDatagrok`,
  `specTestOptions` from
  `public/packages/UsageAnalysis/files/TestTrack/spec-login.ts` —
  used by sibling `radar-spec.ts`.
- **Dataset metadata** (carried over from the original .md trailing
  block): order 28; primary dataset `System:DemoFiles/demog.csv`.
  Step 1 also uses `System:DemoFiles/geo/earthquakes.csv` (declared
  in original Step 1 body but not in the dataset metadata).

---
{
  "order": 28,
  "datasets": ["System:DemoFiles/demog.csv"]
}
