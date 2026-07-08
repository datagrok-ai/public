---
feature: peptides
sub_features_covered:
  - peptides.viewers.sar-base.export-mutation-cliffs
  - peptides.viewers.sar-base
  - peptides.viewers.monomer-position
  - peptides.workflow.start-analysis
target_layer: playwright
coverage_type: smoke
produced_from: atlas-driven
realized_as:
  - export-mutation-cliffs-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:18:30Z
    spec_runs:
      - spec: export-mutation-cliffs-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 123
        failure_keys: []
---

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The `Sequence Variability Map` (MonomerPosition) viewer is added to the TableView alongside the other configured viewers.

## Scenarios

### Scenario 1 — Export Mutation Cliffs opens a new TableView with the documented column shape

Atlas anchor: `critical_paths[export-mutation-cliffs-from-sar-viewer]` (p2) + `interactions[sar-export-mutation-cliffs-to-table]` (smoke). Exercises the context-menu Export surface on `SARViewer` and confirms the exported TableView has the columns produced by `SARViewer.exportMutationCliffs()`.

1. Right-click on the `Sequence Variability Map` viewer body to open the viewer's context menu.
2. Hover the **Export** submenu and click **Export Mutation Cliffs**.
3. In the column-picker dialog, accept the default selection (no extra paired columns) and click **OK**.
4. A new `TableView` is opened (the exported mutation-cliffs table) and becomes the active view.
5. Confirm the new TableView's grid contains the documented columns: `Seq 1`, `Seq 2`, `Mutation` (with semtype `MacromoleculeDifference`, value formatted as `seq1#seq2`), `Seq 1` activity, `Seq 2` activity, `Delta`.
6. Confirm the grid is non-empty (at least one mutation-cliff pair row is present given the peptides.csv demo dataset).

### Scenario 2 — Exported table includes paired columns when extra columns are selected

Same context-menu entry path, but the column-picker selects one additional column to verify the paired-column emission contract.

1. Return to the original peptides TableView (close or switch away from the exported TableView).
2. Right-click on the `Sequence Variability Map` viewer body and choose **Export | Export Mutation Cliffs**.
3. In the column-picker dialog, select one extra non-default column (the first available non-activity column listed) and click **OK**.
4. The new exported `TableView` opens and contains the documented base columns (`Seq 1`, `Seq 2`, `Mutation`, `Seq 1` activity, `Seq 2` activity, `Delta`) plus a paired pair of columns for the extra column chosen (one column carrying the value taken from Seq 1, the other carrying the value taken from Seq 2).

## Notes

- **Atlas provenance.** Scenario 1 derives from `critical_paths[export-mutation-cliffs-from-sar-viewer]` (atlas `peptides.yaml`); the atlas entry's `derived_from:` cites `public/packages/Peptides/src/viewers/sar-viewer.ts#L560`. Scenario 2 derives from `interactions[sar-export-mutation-cliffs-to-table]` (atlas `peptides.yaml`).
- **target_layer rationale.** Multi-step UI flow — viewer context menu, column-picker dialog OK roundtrip, and assertion against the resulting TableView grid shape — requires DOM/UI driving. `apitest` is not viable because the entry point is a context-menu command bound to a Viewer instance and the exported artifact is a new `TableView` (UI surface, not a JS-API return value).
- **Coverage rationale.** `coverage_type: smoke` matches atlas `interactions[sar-export-mutation-cliffs-to-table].coverage_type` verbatim per the canonical rule (atlas value is authoritative when scenario maps onto an atlas interaction/edge entry). The atlas `critical_paths[export-mutation-cliffs-from-sar-viewer].priority: p2` maps via the p0→smoke / p1→regression / p2-p3→author's-discretion heuristic; the dominant atlas declaration here is the `interactions[]` entry's `smoke`, taken as canonical.
- **Coverage contribution.** This scenario retires the largest atlas-anchored uncovered sub_feature in the SAR Export family: `peptides.viewers.sar-base.export-mutation-cliffs` (`viewers.sar-base.export-*` per the chain's Gate F gaps[type: sub-feature-coverage-gap] uncovered family enumeration; see `scenario-chains/peptides.yaml :: gate_f_verdict.gaps[2].uncovered_sub_feature_families`). Sister scenario for the invariant-map export surface (`peptides.viewers.sar-base.export-invariant-map` via atlas `critical_paths[export-invariant-map-from-sar-viewer]` p2 + `interactions[sar-export-invariant-map-to-table]` smoke) remains an open Test Designer authoring item.
- **Setup composition.** Step 2 of Setup launches SAR via the top-menu entry path (`peptides.workflow.sar-dialog` → `peptides.workflow.analyze-ui` → `peptides.workflow.start-analysis`). The context-panel `Launch SAR` button entry path is owned by `sar.md`; this scenario picks the top-menu path so it does not duplicate `sar.md`'s entry-flow coverage and instead focuses the body on the export surface.
- **Related bugs.** Atlas `known_issues[]` carries no curated bug intersecting `peptides.viewers.sar-base.export-mutation-cliffs` or `peptides.viewers.sar-base` as of atlas revision 3; `related_bugs:` is empty by intent.

## Original trailing metadata

```json
{
  "order": 5,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
