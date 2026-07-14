---
feature: peptides
target_layer: playwright
coverage_type: smoke
priority: p2
realizes: [export-invariant-map-from-sar-viewer]
produced_from: atlas-driven
realized_as:
  - export-invariant-map-spec.ts
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
    timestamp: 2026-05-30T21:37:00Z
    spec_runs:
      - spec: export-invariant-map-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 39
        failure_keys: []
---

# Peptides — Export Invariant Map from the SAR viewer

After running SAR analysis, the Sequence Variability Map and Most Potent Residues viewers both offer an **Export Invariant Map** context-menu action that opens a new table of monomer-by-position counts. This scenario checks that the exported table has the right shape (a `Monomer` column plus one column per sequence position, with integer counts) and that the action works identically from either viewer, since both share the same underlying export code.

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The `Sequence Variability Map` (MonomerPosition) viewer is added to the TableView alongside the other configured viewers.

## Scenarios

### Scenario 1 — Export Invariant Map opens a new table with the documented monomer-by-position count shape

Checks the context-menu Export action on the Sequence Variability Map viewer and confirms the exported table has one `Monomer` column plus one column per sequence position, where each cell counts how many sequences carry that monomer at that position.

1. Right-click on the `Sequence Variability Map` viewer body to open the viewer's context menu.
2. Hover the **Export** submenu and click **Export Invariant Map**.
3. A new `TableView` is opened (the exported invariant-map table) and becomes the active view.
4. Confirm the new TableView's grid contains a leading `Monomer` column (the monomer identifier per row).
5. Confirm the new TableView's grid contains one column per sequence position (the number of position columns matches the number of per-position columns the SAR analysis produced on the original TableView).
6. Confirm at least one row is present (the demo dataset's monomer set is non-empty) and that position-column cell values are integers (sequence counts, including zero where no sequence carries the given monomer at that position).
7. Confirm the leading `Monomer` column carries a value visibly recognizable as a monomer identifier (single-letter amino acid or a HELM-shape monomer label, consistent with how monomers render in the original Peptides grid's WebLogo headers).

### Scenario 2 — Export Invariant Map from Most Potent Residues produces the same shape

The Export Invariant Map command is shared by both SAR viewers (Sequence Variability Map and Most Potent Residues), so it must behave identically from either entry point. This scenario repeats the export from Most Potent Residues to confirm that shared behavior empirically, not just by code inspection.

1. Return to the original peptides TableView (close or switch away from any previously exported TableView). The `Most Potent Residues` viewer is present in the TableView (added by the default SAR launch in Setup).
2. Right-click on the `Most Potent Residues` viewer body to open the viewer's context menu.
3. Hover the **Export** submenu and click **Export Invariant Map**.
4. A new `TableView` is opened (the exported invariant-map table).
5. Confirm the exported TableView's grid has the same shape as Scenario 1: a leading `Monomer` column + one column per sequence position; cell values are integer sequence-counts; ≥1 row present.
6. Confirm no error balloon or `Cannot read ... on null` console error is raised during the Export action (the SARViewer-base Export contract must hold regardless of which subclass is the entry point).

## Notes

- **Setup composition.** SAR is launched via the top-menu entry path, mirroring sibling `export-mutation-cliffs.md`. The context-panel `Launch SAR` button entry path is owned by `sar.md`, so this scenario avoids duplicating that coverage.
- **No related bugs.** No curated bug currently anchors the Export Invariant Map surface.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Table View`; `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Sequence Variability Map`.

## Original trailing metadata

```json
{
  "order": 7,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
