---
feature: bio
sub_features_covered:
  - bio.detector
  - bio.rendering
  - bio.transform.convert-notation
  - bio.transform.convert-notation.action
  - bio.io.fasta-handler
  - bio.io.save-as-fasta
  - bio.analyze.sequence-space.transform
  - bio.api.get-seq-helper
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - GROK-12164
  - GROK-15176
  - GROK-18616
  - GROK-19928
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-lifecycle-macromolecule-column-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T00:14:00Z
    spec_runs:
      - spec: bio-lifecycle-macromolecule-column-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 115
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T05:00:00Z
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: NA
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
---

# Bio — `macromolecule_column` source-class lifecycle

Proactive lifecycle scenario for the `macromolecule_column`
source class per
`scenario-chains/bio.yaml#proactive_lifecycle_specs[0]`. Walks
the column through every non-agnostic
`dep_lifecycle_op` declared in atlas
`dep_lifecycle_ops[].affected_source_classes` that touches the
shape: detect on open → convert notation (FASTA → SEPARATOR) →
export as FASTA → re-import FASTA → save project with analysis
output → reopen and verify analysis output survives the
round-trip.

`external_deps: []` per chain `proactive_lifecycle_specs[0]` —
the column is an in-memory `DG.Column` on a `DG.DataFrame`; no
server-side entity-type prerequisite. `env_requirements: []` —
runs against System.AppData test data only.

Distinct from the section's analyze / convert / search runners
(`analyze.md`, `convert.md`, `msa.md`, etc.) in that it covers
the **lifecycle** of the macromolecule column shape (save /
reopen, multi-op chaining) rather than a single-op golden path.
Mitigates the cross-cutting `bio.x.detector-renderer-after-convert`
(GROK-12164), `bio.x.bio-to-chem-via-atomic-level` (GROK-15176
adjacent), `bio.x.entry-path-detector-sync` (GROK-18616), and
`bio.x.bio-analysis-in-datasync-projects` (GROK-19928) atlas
interactions at the lifecycle (save-and-reopen) granularity.

## Setup

- Authenticate to Datagrok as the test user (Playwright session).
- Test dataset: `System.AppData/Bio/tests/filter_FASTA.csv`
  (canonical FASTA Macromolecule column; atlas
  `source_classes[macromolecule_column].examples[0]`).
- Conversion target dataset: `System.AppData/Bio/tests/filter_HELM.csv`
  (HELM Macromolecule column — used as Step 2's source notation
  to exercise the HELM→SEPARATOR convert branch that hit
  GROK-12164).
- Project name: `bio-lifecycle-macromolecule-${Date.now()}`
  (unique per run to keep parallel-run safe; cleaned in Step 6).
- Service-surface dependency: `await
  grok.functions.call('Bio:getSeqHelper')` is used in Step 2 and
  Step 4 as the supported entry-point for column manipulation
  (atlas `bio.api.get-seq-helper`).
- Cleanup: Step 6 deletes the project and removes the exported
  FASTA file.

## Scenarios

### Scenario 1: Detect on open + convert-notation round trip

Steps:
1. Open `System.AppData/Bio/tests/filter_HELM.csv` via
   `grok.data.loadTable` (or Browse > Files double-click). The
   Macromolecule detector (atlas `bio.detector`) classifies the
   sequence column synchronously on open; the `units` /
   `aligned` / `alphabet` / `separator` column tags are set
   before any analyze dialog opens.
2. Verify the detector outcome:
   - `column.semType === 'Macromolecule'`.
   - `column.getTag('units') === 'helm'`.
   - The HELM cell renderer (atlas `bio.rendering` →
     `bio.rendering.biln` / `bio.rendering.custom` per detected
     unit) paints sequence cells (no `[object Object]` fallback,
     no blank canvas).
3. Right-click the Macromolecule column header > Convert ... >
   pick **SEPARATOR** (atlas
   `bio.transform.convert-notation.action`,
   `bio.transform.convert-notation`). Use the dash `-` as
   separator — the same notation pair that surfaced GROK-12164.
4. Verify the convert outcome:
   - A new SEPARATOR-units Macromolecule column appears.
   - `column.getTag('units') === 'separator'`.
   - The separator cell renderer (atlas `bio.rendering.separator`)
     dispatches and paints correctly divided monomers (no
     monomer-collapse, no missing dividers).

Expected:
- Detector tags are stable on open (no late mutation; reading
  `units` immediately after load returns the final value).
- Convert mutates the column shape and the renderer-dispatch
  follows the new tags consistently — the
  `bio.x.detector-renderer-after-convert` cross-feature contract
  (GROK-12164) is satisfied at the lifecycle layer.

### Scenario 2: Import FASTA → Export FASTA → re-import round trip

Steps:
1. Open `System.AppData/Bio/tests/filter_FASTA.csv` (canonical
   FASTA Macromolecule sample). The Bio FASTA file handler
   (atlas `bio.io.fasta-handler`) ingests; the detector
   classifies synchronously.
2. Right-click the Macromolecule column > Save As FASTA... (or
   File > Export > As FASTA...; atlas `bio.io.save-as-fasta`).
   Save to a temp path (e.g.
   `System.AppData/UsageAnalysis/temp/lifecycle-${Date.now()}.fasta`).
3. Re-open the exported FASTA file via `grok.data.loadTable`
   against the temp path. The FASTA handler runs again on a file
   that this scenario itself produced — round-trippable contract
   (atlas `bio.cp.import-fasta-export-fasta`,
   `bio.cp.fasta-import-via-multiple-entry-paths`).
4. Verify the imported table:
   - Row count matches the original FASTA's sequence count
     (modulo any silent dedupe — assert exact equality).
   - The first re-imported sequence, normalized, equals the
     first original sequence (string equality after
     `getSeqHelper`-aware comparison).
   - The Macromolecule renderer dispatches the FASTA cell
     renderer (atlas `bio.rendering.fasta`).

Expected:
- The exported FASTA is round-trippable: re-importing reproduces
  the original sequence content (no monomer-name drift, no
  separator artifacts, no truncation).
- The entry-path detector-sync invariant holds for the
  programmatic load path (`grok.data.loadTable`) per atlas
  `bio.x.entry-path-detector-sync` (GROK-18616) — the dialog
  side is exercised in Scenario 3 via Sequence Space.

### Scenario 3: Save project with analysis + reopen restores analysis output

Steps:
1. With the FASTA table from Scenario 2 active (or re-open
   `filter_FASTA.csv` if Scenario 2 was not run), trigger
   **Bio | Analyze | Sequence Space...** (atlas
   `bio.analyze.sequence-space.top-menu`,
   `bio.analyze.sequence-space.editor`). Run with default
   parameters. Wait for the embedding compute to finish; the
   ScatterPlot embedding viewer docks (atlas
   `bio.analyze.sequence-space.transform`).
2. Verify the embedding output is present:
   - The two new embedding columns (`Embed_X`, `Embed_Y` or
     equivalent per atlas
     `bio.analyze.sequence-space.transform`) are added to the
     DataFrame.
   - The ScatterPlot viewer is in the active view's viewer
     list.
3. Save the project: Ribbon SAVE button (NOT Ctrl+S, per the
   `feedback_no_ctrlS_for_layouts` policy); project name from
   Setup; **Data Sync** toggle ON. Cancel the auto-share dialog
   if it appears. Verify `POST /projects` succeeds and
   `grok.dapi.projects.find(<id>)` returns the saved project.
4. Close the project / reopen via
   `grok.dapi.projects.find(<id>)` → `project.open()` (or the
   equivalent JS API path; UI driving for project reopen is
   delegated to UI-smoke scenarios elsewhere — this scenario
   asserts the persistence-side outcome).
5. Verify the reopened project state:
   - The Macromolecule column is back with the same `units` /
     `semType` tags (detector state survives the save-reopen
     round trip).
   - The Sequence Space embedding columns are present.
   - The ScatterPlot viewer (or its persisted layout) restores;
     `isRowDrawable` style invariants hold (no NaN / missing
     coordinate rows that silently drop from the plot).

Expected:
- Project save + reopen with Data Sync ON preserves the Bio
  analysis output shape — the
  `bio.x.bio-analysis-in-datasync-projects` cross-feature
  contract (GROK-19928) is satisfied at the lifecycle layer.
- No silent persistence drop on reopen (embedding columns
  intact, viewer restored).

### Scenario 4: Cleanup

Steps:
1. Delete the project:
   `await grok.dapi.projects.delete(project)`.
2. Delete the temp FASTA file produced by Scenario 2 (best-
   effort; ignore not-found errors).
3. Close any docked viewers / dialogs left over.

Expected:
- Cleanup runs in `tearDownAll` / `finally` regardless of
  earlier failures (no test-state leak across runs per the
  testing-rules in `core/docs/platform/TESTING.md`).

## Notes

- Origin: chain rev 6
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[0]`
  (`source_class: macromolecule_column`, `spec_target:
  bio-lifecycle-macromolecule-column-spec.ts`). Lands the
  intent-layer `.md` realization that
  `F-PROACTIVE-COVERAGE-01` (SCOPE_REDUCTION, cycle
  2026-05-31-bio-migrate-01) cited as the gap. `.ts`
  realization (`spec_target`) is downstream Automator / Phase
  C, not gated by F.
- atlas entry derived from
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[0]`
  (chain author, cycle 2026-05-30-bio-chain-bootstrap-01).
- `dep_lifecycle_ops` exercised (per
  `proactive_lifecycle_specs[0].bundled_ops` minus
  `save_monomer_library` / `load_monomer_library` /
  `align_sequences` which do not affect macromolecule_column):
  `detect_macromolecule_on_open` (Scenario 1.1, 2.1),
  `convert_notation` (Scenario 1.3),
  `import_fasta_file` (Scenario 2.1, 2.3),
  `export_as_fasta` (Scenario 2.2),
  `save_project_with_analysis` (Scenario 3.3).
- target_layer rationale: `playwright` — the lifecycle spans
  top-menu dispatch (Bio | Analyze | Sequence Space), context-
  menu Convert action, file-handler ingest, Save Project
  dialog, and viewer-docking verification. `apitest` cannot
  exercise the ribbon path nor the Save Project UI; JS API
  substitutes are used for the persistence-side assertions
  (Step 3.4, 4.1) per the same pattern as sibling
  `projects-lifecycle-*.md` scenarios — UI driving stays on
  the dispatch points where the assertable surface lives.
- Sibling tests considered:
  - `public/packages/Bio/src/tests/fasta-export-tests.ts`
    covers `saveAsFastaDo` round-trip at the API layer (atlas
    `bio.io.save-as-fasta`); this scenario adds the missing
    UI-layer + project-persistence lifecycle layer.
  - `public/packages/Bio/src/tests/sequence-space-tests.ts`
    covers the embedding compute API (atlas
    `bio.analyze.sequence-space.transform`); this scenario
    adds the save-and-reopen persistence layer that the
    apitest does not exercise.
  - `public/packages/Bio/src/tests/detectors-tests.ts` /
    `convert-tests.ts` cover detector + convert at the unit
    layer; this scenario adds the cross-op chaining layer.
- This scenario covers 8 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2 — well above;
  `F-STRUCT-INTERACTION-01` floor: 3 — satisfied per scenario).
  Scenario cardinality (per `## Scenarios` section): 4 (one
  cleanup + three substantive lifecycle scenarios) — meets the
  ≥ 2 scenarios floor.
- Manual-only subset: none of the eight covered sub_features
  appear in atlas `manual_only[]` (verified against atlas rev 3
  `manual_only[]` list: `bio.viewers.web-logo`,
  `bio.viewers.vd-regions`, `bio.rendering.column-header`,
  `bio.rendering.macromolecule-difference`, the demo entries,
  `bio.panels.{structure-3d, atomic-level, tooltip}`). The
  scenario does not claim coverage over
  `bio.rendering.macromolecule-difference` even though the
  convert-notation invariant is adjacent to GROK-16596 reset-
  state behavior — the diff-cell renderer is manual-only and
  not covered here.
- `coverage_type: regression` per STEP E heuristic: this is
  general coverage of the macromolecule_column lifecycle shape
  (not a single critical_path golden path → not smoke; not a
  boundary value → not edge; not stress/latency-sensitive →
  not perf). Atlas does not register a specific edge_cases[]
  entry that maps onto this multi-op lifecycle; STEP E
  fallback to `regression` applies.
- Deferrals: none mandatory. Step 3.4's UI reopen path uses JS
  API by design (consistent with sibling
  `projects-lifecycle-*.md` scenarios); driving the full UI
  reopen (Browse > Dashboards > tile click) is a UI-smoke
  responsibility elsewhere in the section.
- Related-bug context: `related_bugs: [GROK-12164, GROK-15176,
  GROK-18616, GROK-19928]` per chain
  `proactive_lifecycle_specs[0].bugs_reinforcing` (GROK-15176
  reinforced via the atomic-level / molfile-validity branch
  that an atlas `bio.x.bio-to-chem-via-atomic-level` lifecycle
  cell would reach; this scenario reinforces it at the convert
  / round-trip persistence layer without claiming the cross-
  package PubChem contract — that is the bug-focused spec
  `bio-grok-15176-spec.ts` per chain
  `bug_focused_candidates`).

---
{
  "order": 13
}
