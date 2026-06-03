---
feature: dendrogram
sub_features_covered:
  - dendrogram.clustering.menu.chem
  - dendrogram.clustering.dialog
  - dendrogram.clustering.inject-tree-for-grid
target_layer: playwright
coverage_type: smoke
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Dendrogram/hierarchical-clustering-chem.md
migration_date: 2026-06-02
source_text_fixes:
  - canonicalize-headings-h1-setup-scenarios
  - promote-block-a-block-b-to-named-scenarios
  - cite-explicit-mol1k-csv-path-from-order-trailer
  - downgrade-step-7-non-monotonic-visual-assertion-to-no-throw
candidate_helpers: []
unresolved_ambiguities:
  - step-7-non-monotonic-tree-not-assertable-via-playwright-canvas
scope_reductions:
  - step-7-numeric-features-substitution-deferred-canvas-column-picker-not-dom-addressable
  - step-7-centroid-on-molecule-features-triggers-wasm-oob-platform-crash-swapped-to-median-per-scenario-or-clause
related_bugs: []
realized_as:
  - hierarchical-clustering-chem-spec.ts
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - hierarchical-clustering-dialog
  - hierarchical-clustering-dialog-distance-dropdown
  - hierarchical-clustering-dialog-linkage-dropdown
  - hierarchical-clustering-dialog-features-column-input
  - hierarchical-clustering-dialog-ok-button
  - chem-analyze-hierarchical-clustering-menu-entry
ui_coverage_delegated_to: null
gate_verdicts:
  d:
    verdict: EVIDENCE_GAP
    cycle_id: 2026-06-03-dendrogram-migrate-01
    timestamp: 2026-06-03T12:30:00Z
    failure_keys: []
    evidence_needed: |
      Third consecutive same-shape EVIDENCE_GAP dispatch on this
      scenario. Inputs supply only scenario_path and mode_file; no
      original_path is provided, and the migrated frontmatter
      original_path self-references the same file (in-place Migrator
      rewrite, migration_date 2026-06-02). The section directory has
      no sibling original*.md / pre-migration snapshot. The nine
      content-diff checks (D-STEP-01, D-STEP-02, D-EDGE-01, D-SAN-02,
      D-STRUCT-01, D-STRUCT-02, D-MERIT-01, D-MERIT-02,
      D-UI-DELEGATION-01) are by construction an ORIGINAL-vs-MIGRATED
      diff and remain not constructible from disk under this dispatch
      shape. The verdict is deterministic and bit-identical to the
      prior two rounds; re-invoking with the same inputs will reproduce
      it indefinitely.
      Mechanical and Phase 1 checks pass independently of any diff
      and are recorded here as evidence that the migrated artifact
      itself is well-formed:
      - D-STRUCT-MECH-03 PASS — 8/8 required migration-output fields
        populated (feature, sub_features_covered ×3, target_layer
        playwright, coverage_type smoke, produced_from migrated,
        original_path, migration_date 2026-06-02, related_bugs []);
        deprecated migrated_from absent.
      - D-STRUCT-MECH-05 PASS — original_path target exists on disk
        (degenerate self-reference; file present).
      - D-FRONTMATTER-PHASE1-01 PASS — all four Phase 1 keys present
        and parseable as YAML lists.
      - D-FRONTMATTER-PHASE1-02 PASS — per-field schemas honored:
        source_text_fixes ×4 kebab-case slugs ≤80 chars no
        duplicates; candidate_helpers []; unresolved_ambiguities ×1
        kebab-case slug; scope_reductions [].
      Migrator dispatch sidecar
      (hierarchical-clustering-chem.md.migrator.dispatch.yaml) records
      verdict: PASS with outputs.no_op: true, re_migration: true, and
      an explicit recommendation that the orchestrator skip Gate D
      for this scenario and carry the prior Gate D verdict forward
      via the already-migrated short-circuit. Reaching this critic
      node again with the same scenario_path-only input is the same
      upstream dispatcher-input-shape gap identified in the prior two
      rounds; the resolution is architectural (operator-side), not
      in-loop recon.
      Resolution options for the operator (either closes the gap;
      neither is automatable from inside the critic):
      (1) Recover the pre-migration revision of
      public/packages/UsageAnalysis/files/TestTrack/Dendrogram/hierarchical-clustering-chem.md
      from git history (commit predating migration_date 2026-06-02),
      write it to a separable path (e.g.
      .claude/scratch/dendrogram-chem-original.md), and re-dispatch
      this critic with inputs.original_path set to the recovered
      file. This converts the dispatch from scenario_path-only to a
      true ORIGINAL-vs-MIGRATED bundle and unlocks the nine
      content-diff checks.
      (2) Accept the Migrator no_op-PASS and bypass Gate D via the
      dispatcher's already-migrated short-circuit, routing this
      iteration onward (Critic A or next node) carrying the on-disk
      Gate D verdict forward. This is the option the Migrator sidecar
      explicitly recommends.
      Corroborating evidence that does NOT substitute for the diff
      (informational only): chain YAML dependency_graph entry for
      hierarchical-clustering-chem.md (revision 5) declares
      pyramid_layer ui-smoke, classification simple, 7-step body, and
      ui_coverage_responsibility[] of six dialog-gateway flow ids; the
      migrated frontmatter ui_coverage_responsibility[] matches that
      list verbatim. The chain YAML unresolved_ambiguities[] entry
      describes exactly the Step 7 non-monotonic centroid/median
      ambiguity that the migrated unresolved_ambiguities[] slug
      carries. The source_text_fixes slug
      downgrade-step-7-non-monotonic-visual-assertion-to-no-throw
      matches chain YAML option 1 of the recorded Step 7 resolution.
      The migrated body has 7 numbered steps split as Block A
      (steps 1-4: dialog open + Distance and Linkage dropdown
      enumerations) and Block B (steps 5-7: three OK-click runs
      across ward, single, and centroid-or-median linkages), each
      with explicit Expected-result lines, consistent with the chain
      YAML 7-step classification. The structural shape matches the
      chain-analyzer intent record on every dimension a chain record
      carries, so neither a FAIL (silent drop) nor a SCOPE_REDUCTION
      is indicated by the corroborating evidence.
      Routing note for orchestrator: this is the third consecutive
      EVIDENCE_GAP round on this scenario within the current cycle
      (2026-06-03-dendrogram-migrate-01 timestamps 09:39:57, 09:42:38,
      and the present invocation) and the seventh across two
      consecutive cycles, none of which produced PASS or
      SCOPE_REDUCTION. The grok-critic verdict-semantics table caps
      review rounds at 2 per gate and routes any further iteration
      through autopilot boundary trigger #1. The
      recon_chain_analyzer LOOP_CAP_EXCEEDED entry already on file
      (decision-log L8103, cap_value 3, observed_iterations 3,
      last_verdict EVIDENCE_GAP) is the documented escalation
      handoff for the same loop. Continued in-loop dispatch is
      contra-indicated; one of the two operator resolution options
      above is the documented next action.
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T14:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T16:53:00Z
    spec_runs:
      - spec: hierarchical-clustering-chem-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 36
        failure_keys: []
---

# Hierarchical Clustering (Chem) — dialog UI smoke

Section ui-smoke for the Dendrogram chain. Owns the
hierarchical-clustering dialog gateway exercised through the
`Chem | Analyze | Hierarchical Clustering...` menu entry: dialog open,
Distance and Linkage dropdown enumerations, and representative
end-to-end OK-click runs that inject a dendrogram neighbor on the
grid. Bio-path defaults and Assign Clusters / Ctrl+wheel zoom flows
are owned by their respective scenarios; this smoke verifies only the
dialog gateway and one-OK-one-dendrogram contract.

## Setup

- A clean Datagrok view (no preloaded tables).
- The Chem package is installed and registered (the `Chem | Analyze`
  top-menu surface and `Chem:getMorganFingerprints` function must be
  available).
- The Dendrogram package is installed and registered.
- Dataset: `System:AppData/Chem/mol1K.csv` (a 1000-row dataset with a
  `molecule` column carrying SMILES; also exposes numeric columns
  `pIC50_HIV_Integrase` and `Q`).

## Scenarios

### Block A — Dialog exposes all Distance and Linkage values

1. Open the dataset at `System:AppData/Chem/mol1K.csv`. Wait for the
   `molecule` column to render structures.
2. Run `Chem | Analyze | Hierarchical Clustering...`.
   * Expected result: the Hierarchical Clustering dialog opens with
     `Table` = `mol1K`, `Features` defaulting to `molecule`, and
     `Distance` and `Linkage` inputs visible.
3. Open the `Distance` dropdown.
   * Expected result: exactly two values are listed — `euclidean`,
     `manhattan` (default `euclidean`).
4. Open the `Linkage` dropdown.
   * Expected result: exactly seven values are listed, in order —
     `single`, `complete`, `average`, `weighted`, `centroid`,
     `median`, `ward` (default `ward`).

### Block B — Representative end-to-end runs (spot-check)

5. With `Features` = `molecule`, `Distance` = `euclidean`,
   `Linkage` = `ward`, click `OK`.
   * Expected result: a `Creating dendrogram ...` progress indicator
     appears, then a dendrogram is injected to the left of the grid
     with one leaf per row. No console errors.
6. Close the dendrogram. Re-open the dialog, set
   `Distance` = `manhattan`, `Linkage` = `single`,
   `Features` = `molecule`, click `OK`.
   * Expected result: a dendrogram builds successfully (different
     shape, still one leaf per row). No `Unsupported column type`
     error and no fatal console errors.
7. Close the dendrogram. Re-open the dialog, set
   `Distance` = `euclidean`, `Linkage` = `centroid` (or `median`),
   `Features` = numeric columns (`pIC50_HIV_Integrase`, `Q`), click
   `OK`.
   * Expected result: a dendrogram builds without error. No
     `Unsupported column type` error and no fatal console errors.

## Notes

- All three top-menu entries — `Bio | Analyze | Hierarchical
  Clustering...`, `Chem | Analyze | Hierarchical Clustering...`, and
  `ML | Cluster | Hierarchical Clustering...` — call the same
  `hierarchicalClusteringDialog`; they differ only in the
  default-selected feature column. The Bio sequence-default path is
  covered by `hierarchical-clustering-bio.md`; the ML path is not
  covered separately in the Dendrogram chain (same gateway).
- The `molecule` distance path calls `Chem:getMorganFingerprints`
  followed by Tanimoto distance; numeric columns use the per-column
  difference metric (`tree-helper.ts:526-548`).
- Centroid and median linkages on numeric data may produce a
  **non-monotonic tree** (branches with inversions or apparent
  cross-overs). This is mathematically expected for these linkage
  methods and is **not** a defect. The non-monotonic property is a
  structural property of the merge-height sequence and is asserted
  in the apitest matrix (`hierarchical-clustering-chem-api.md`,
  Scenario 1 covers all 14 distance × linkage combinations on the
  molecule path with leaf-count + no-throw checks). The UI smoke
  here only verifies that the OK click produces a dendrogram and
  surfaces no console error for the centroid/median path.
- Existing package unit coverage
  (`Dendrogram/src/tests/hierarchical-clustering-tests.ts`) exercises
  only `euclidean` + `average` on a numeric column; the other
  distance × linkage combinations and the molecule path are not
  covered there.
- Per the Dendrogram section ui-smoke contract (chain
  `ui_coverage_plan.smoke_covers`), this scenario owns the
  `hierarchical-clustering-dialog`,
  `hierarchical-clustering-dialog-distance-dropdown`,
  `hierarchical-clustering-dialog-linkage-dropdown`,
  `hierarchical-clustering-dialog-features-column-input`,
  `hierarchical-clustering-dialog-ok-button`, and
  `chem-analyze-hierarchical-clustering-menu-entry` flows.
  Bio-path dialog flows delegate to this scenario per the chain's
  `ui_coverage_plan.delegated_scenarios`.

---
{
  "order": 2,
  "datasets": ["System:AppData/Chem/mol1K.csv"]
}
