---
feature: peptides
sub_features_covered:
  - peptides.viewers.sar-base.export-invariant-map
  - peptides.viewers.sar-base
  - peptides.viewers.monomer-position
  - peptides.workflow.start-analysis
target_layer: playwright
coverage_type: smoke
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

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The `Sequence Variability Map` (MonomerPosition) viewer is added to the TableView alongside the other configured viewers.

## Scenarios

### Scenario 1 — Export Invariant Map opens a new TableView with the documented monomer-by-position count shape

Atlas anchor: `critical_paths[export-invariant-map-from-sar-viewer]` (p2) + `interactions[sar-export-invariant-map-to-table]` (smoke). Exercises the context-menu Export surface on `SARViewer` and confirms the exported TableView has the columns produced by `SARViewer.exportInvariantMap()` (one `Monomer` column + one column per sequence position; cell value = count of sequences carrying that monomer at that position).

1. Right-click on the `Sequence Variability Map` viewer body to open the viewer's context menu.
2. Hover the **Export** submenu and click **Export Invariant Map**.
3. A new `TableView` is opened (the exported invariant-map table) and becomes the active view.
4. Confirm the new TableView's grid contains a leading `Monomer` column (the monomer identifier per row).
5. Confirm the new TableView's grid contains one column per sequence position (the number of position columns matches the number of per-position columns the SAR analysis produced on the original TableView).
6. Confirm at least one row is present (the demo dataset's monomer set is non-empty) and that position-column cell values are integers (sequence counts, including zero where no sequence carries the given monomer at that position).
7. Confirm the leading `Monomer` column carries a value visibly recognizable as a monomer identifier (single-letter amino acid or a HELM-shape monomer label, consistent with how monomers render in the original Peptides grid's WebLogo headers).

### Scenario 2 — Export Invariant Map from Most Potent Residues viewer produces the same shape (SARViewer-base contract)

Atlas anchor: `peptides.viewers.sar-base` is the abstract base for both `peptides.viewers.monomer-position` and `peptides.viewers.most-potent-residues`; the Export Invariant Map command is registered on the base class and so must be available from both viewer subclasses with the same exported-table shape. This scenario exercises the contract from the `Most Potent Residues` entry point so the SARViewer-base inheritance of the Export action is empirically covered, not just inferred.

1. Return to the original peptides TableView (close or switch away from any previously exported TableView). The `Most Potent Residues` viewer is present in the TableView (added by the default SAR launch in Setup).
2. Right-click on the `Most Potent Residues` viewer body to open the viewer's context menu.
3. Hover the **Export** submenu and click **Export Invariant Map**.
4. A new `TableView` is opened (the exported invariant-map table).
5. Confirm the exported TableView's grid has the same shape as Scenario 1: a leading `Monomer` column + one column per sequence position; cell values are integer sequence-counts; ≥1 row present.
6. Confirm no error balloon or `Cannot read ... on null` console error is raised during the Export action (the SARViewer-base Export contract must hold regardless of which subclass is the entry point).

## Notes

- **Atlas provenance.** Scenario 1 derives from `critical_paths[export-invariant-map-from-sar-viewer]` (atlas `peptides.yaml`); the atlas entry's `derived_from:` cites `public/packages/Peptides/src/viewers/sar-viewer.ts#L674`. Scenario 2 derives from the same `critical_paths[]` entry + `interactions[sar-export-invariant-map-to-table]` (atlas `peptides.yaml`); the SARViewer-base inheritance contract is rooted at `public/packages/Peptides/src/viewers/sar-viewer.ts#L104` (`peptides.viewers.sar-base` atlas source line).
- **target_layer rationale.** Multi-step UI flow — viewer context menu, Export submenu navigation, new-TableView assertion against grid shape (column count, column names, cell types) — requires DOM/UI driving. `apitest` is not viable because the entry point is a context-menu command bound to a Viewer instance and the exported artifact is a new `TableView` (UI surface, not a JS-API return value). Same layer choice as the realized sister scenario `export-mutation-cliffs.md`.
- **coverage_type rationale.** Atlas `interactions[sar-export-invariant-map-to-table].coverage_type: smoke` is canonical per the rule "When the scenario being authored maps onto an atlas `interactions[]` entry, the frontmatter `coverage_type:` MUST match that entry's `coverage_type:` verbatim. Do NOT re-derive." Scenario carries `coverage_type: smoke` to match the atlas entry. The atlas `critical_paths[export-invariant-map-from-sar-viewer].priority: p2` (severity axis) is consistent with the smoke test-kind via the p2→author's-discretion mapping; the dominant atlas declaration here is the `interactions[]` entry's `smoke`.
- **Coverage contribution.** This scenario retires `peptides.viewers.sar-base.export-invariant-map` (previously uncovered, on the `viewers.sar-base.export-invariant-map` arm of the Gate F gaps[type: sub-feature-coverage-gap] uncovered-family enumeration; see `scenario-chains/peptides.yaml :: gate_f_verdict.gaps[1].uncovered_sub_feature_families[viewers.{position-statistics, sar-base.export-invariant-map}]`). Sister-of-realized scenario for the SAR Export family (the mutation-cliffs export arm was realized as `export-mutation-cliffs.md`). The other three `sub_features_covered[]` entries — `peptides.viewers.sar-base`, `peptides.viewers.monomer-position`, `peptides.workflow.start-analysis` — are intentional re-coverage to satisfy `F-STRUCT-INTERACTION-01` (≥3 sub_features per scenario) on the same SAR-Export shape.
- **Related bugs.** Atlas `known_issues[]` carries no curated bug intersecting `peptides.viewers.sar-base.export-invariant-map` or `peptides.viewers.sar-base` as of atlas revision 3 (parallel to the realized `export-mutation-cliffs.md`'s empty `related_bugs:`); `related_bugs:` is empty by intent.
- **Setup composition.** Step 2 of Setup launches SAR via the top-menu entry path (`peptides.workflow.sar-dialog` → `peptides.workflow.analyze-ui` → `peptides.workflow.start-analysis`). The context-panel `Launch SAR` button entry path is owned by `sar.md`; this scenario picks the top-menu path to mirror the sister `export-mutation-cliffs.md` and avoid duplicating `sar.md`'s entry-flow coverage.
- **Cross-scenario shape coverage.** Scenario 1 exercises the Export action from `Sequence Variability Map` (the `MonomerPosition` subclass — the more commonly used entry point); Scenario 2 exercises it from `Most Potent Residues` (the `MostPotentResidues` subclass). Together they cover the SARViewer-base Export contract from both subclass entry points without requiring a separate dedicated scenario per subclass.
- **Deferrals.** None. Scenario steps map to existing Playwright-driveable surfaces (top-menu invocation, viewer context-menu navigation, exported TableView grid-shape assertion).
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Table View` (atlas help_docs rich-object form maps the `Table View` heading to `peptides.workflow.start-analysis`); `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Sequence Variability Map` (maps to `peptides.viewers.monomer-position`).

## Original trailing metadata

```json
{
  "order": 7,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
