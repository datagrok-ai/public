---
feature: peptides
target_layer: playwright
coverage_type: smoke
priority: p2
realizes_atlas: [export-mutation-cliffs-from-sar-viewer]
realizes: []
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

# Peptides — Export Mutation Cliffs from the SAR viewer

After running SAR analysis, the Sequence Variability Map viewer's context menu offers an **Export Mutation Cliffs** action that opens a new table listing pairs of sequences differing by a single mutation, along with their activity delta. This scenario checks the exported table's column shape, and that extra columns chosen in the export dialog are emitted as paired Seq-1/Seq-2 values.

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The `Sequence Variability Map` (MonomerPosition) viewer is added to the TableView alongside the other configured viewers.

## Scenarios

### Scenario 1 — Export Mutation Cliffs opens a new table with the documented column shape

Checks the context-menu Export action on the Sequence Variability Map viewer and confirms the exported table has the expected columns.

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

- **Setup composition.** SAR is launched via the top-menu entry path so this scenario does not duplicate `sar.md`'s context-panel entry-flow coverage and can focus on the export surface.
- **No related bugs.** No curated bug currently anchors the Export Mutation Cliffs surface.

## Original trailing metadata

```json
{
  "order": 5,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
