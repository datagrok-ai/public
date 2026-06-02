---
feature: bio
sub_features_covered:
  - bio.io.fasta-handler
  - bio.io.save-as-fasta
  - bio.detector
  - bio.rendering
  - bio.rendering.fasta
  - bio.api.get-seq-helper
  - bio.lifecycle.init
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - GROK-18616
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-lifecycle-fasta-file-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T07:00:00Z
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
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T06:20:00Z
    failure_keys: []
    scope_reduction_proposal: |-
      Two scenario steps are realized via JS API substitution rather than the UI affordance the scenario step text names; both substitutions have a sound technical rationale documented in the spec body header (lines 88-113) but are NOT surfaced in scenario .md frontmatter scope_reductions[] (currently []). Migrator must apply the following two entries. SR-01 id=save-as-fasta-dialog-jsapi-substitution check=E-TRACE-02 (S3.1) rationale=Scenario Step 3.1 names File | Export | As FASTA... ribbon path / right-click Save As FASTA... column-picker dialog. bio.md grok-browser reference documents no [name=...] selectors for the Save As FASTA column-picker dialog (verified by grep on .claude/skills/grok-browser/references/bio.md). Sibling test public/packages/Bio/src/tests/fasta-export-tests.ts (saveAsFastaTest1) precedent uses SeqHandler.getSplitted as the wrapped primitive, which is the same code path saveAsFastaUI dispatches into post-OK. Spec exercises this primitive (line 488-501) plus dapi.files.writeAsText (line 510) plus Bio:importFasta re-import (line 522) to assert the atlas bio.cp.import-fasta-export-fasta round-trippable contract. verdict_status=SCOPE_REDUCTION. SR-02 id=save-project-ribbon-jsapi-substitution check=E-TRACE-02 (S4.1) rationale=Scenario Step 4.1 names the ribbon SAVE button (Data Sync ON, cancel auto-share dialog). The Save Project ribbon dialog is a platform-wide UI surface not documented in bio.md selectors; UI driving is delegated to platform-side ui-smoke scenarios elsewhere. Spec uses the canonical test-area helper helpers/projects.ts saveAllTablesWithProvenance (mirrors sibling bio-lifecycle-macromolecule-column-spec.ts S3.3-3.5 and Projects/projects-lifecycle-files-spec.ts), preserving the assertable contract (persisted project + reopened detector state per atlas dep_lifecycle_ops[save_project_with_analysis]). verdict_status=SCOPE_REDUCTION. After Migrator applies these SR entries, the spec is otherwise PASS-compliant. Mechanical checks pass: E-STRUCT-MECH-01 spec exists; -02 parses; -03 test() block at line 178; -04 test name 'Bio fasta_file source-class lifecycle' substring-matches scenario heading 'Bio - fasta_file source-class lifecycle'; -05 imports from @playwright/test + sibling test-area relative imports ../spec-login + ../helpers/projects (test-area relative path convention used across 33+ sibling specs); -06 leading /* --- sub_features_covered: [...] --- */ block present at lines 1-10, all 7 ids (bio.io.fasta-handler, bio.io.save-as-fasta, bio.detector, bio.rendering, bio.rendering.fasta, bio.api.get-seq-helper, bio.lifecycle.init) resolve in feature-atlas/bio.yaml. E-SEL-01/02 both [name=...] selectors are class-1 in bio.md ([name=div-Bio] L76/L606; [name=viewer-Grid] standard platform, also at L107 for grid-canvas warm-up). E-SEL-03 no .fill() calls (none used). E-HELP-01 softStep + loginToDatagrok registered in helpers-registry.yaml lines 3927-3938; saveAllTablesWithProvenance/reopenAndAssertProvenance/deleteProjectWithCleanup exist as canonical TestTrack-area helpers at helpers/projects.ts:926,999,1060 used across 33+ section specs (advisory: registry maintenance gap, not Automator-introduced reinvention). E-HELP-02 no reinvention. E-LAYER-01 target_layer=playwright matches paradigm; E-LAYER-COMPLIANCE-01 satisfied by 2 page.locator(...).waitFor calls at lines 240 + 244. E-BOUND-01/02 spec under public/packages/UsageAnalysis/files/TestTrack/bio/, no core/** edits. E-RETRY-IGNORES-GATE-B: retry context active (gate_verdicts.b.verdict=FAIL same-cycle 2026-06-01-bio-migrate-02 with failure_keys [B-RUN-PASS, B-STAB-01]); new spec substantively addresses Validator's diagnosis - corrected {content} -> {fileContent} at 3 call sites (lines 218, 396, 522) per Playwright report error 'Errors calling importFasta: fileContent: Value not defined' verbatim cited in retry dispatch; moved FASTA open into setup phase OUTSIDE softStep (lines 212-239) so failures surface immediately instead of being swallowed; removed broken grok.data.loadTable(System:AppData/...) fallback; added grok.data.detectSemanticTypes(df) at 3 sites (lines 227, 400, 525); evidence-backed test-bug fix triangulated from package.ts:1107-1118 decorator + renderers-test.ts:134 sibling pattern. Not cosmetic, not paradigm-pivot-without-backing, not ignored. E-RETRY-IGNORES-GATE-B does NOT fire.
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T05:40:00Z
    spec_runs:
      - spec: bio-lifecycle-fasta-file-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 50
        failure_keys: []
---

# Bio — `fasta_file` source-class lifecycle

Proactive lifecycle scenario for the `fasta_file` source class
per `scenario-chains/bio.yaml#proactive_lifecycle_specs[3]`.
Walks a FASTA file through every non-agnostic
`dep_lifecycle_op` declared in atlas
`dep_lifecycle_ops[].affected_source_classes` that touches the
shape: drop / open the file via multiple entry paths → the
FASTA handler converts to Macromolecule column → detector
classifies synchronously → export As FASTA round trip → save a
project that includes the imported table → reopen and verify
the FASTA-imported shape survives the round trip.

`external_deps: [FileShare]` per chain `proactive_lifecycle_specs[3]`
— the lifecycle reads/writes under `System:AppData/Bio/samples`
and a temp path for the export. `env_requirements: ["FileShare
connection with per-user permissions"]`.

Mitigates `GROK-18616` (entry-path-class detection-sync gap) by
exercising the multiple-entry-path lifecycle shape from atlas
critical_path `bio.cp.fasta-import-via-multiple-entry-paths`.
Distinct from `bio-lifecycle-macromolecule-column.md` (which
covers a single import path inside its Scenario 2 round trip)
in that this scenario explicitly exercises **multiple** entry
paths within the same run, asserting detector-sync invariants
hold across all of them.

## Setup

- Authenticate to Datagrok as the test user (Playwright session)
  with read/write access to `System:AppData/Bio/samples` and a
  temp directory for exports.
- Source FASTA file: any of `human_genes.fasta` (per atlas
  `source_classes[fasta_file].examples[0]`, the
  `GROK-18474`-tagged sample). If absent, fall back to any
  `*.fasta` file under `System:AppData/Bio/samples`.
- Export target path:
  `System:AppData/UsageAnalysis/temp/lifecycle-fasta-${Date.now()}.fasta`.
- Project name:
  `bio-lifecycle-fasta-file-${Date.now()}`.
- Service-surface dependency:
  `await grok.functions.call('Bio:getSeqHelper')` is used in
  Steps 1 and 2 (atlas `bio.api.get-seq-helper`).
- Cleanup: Step 5 deletes the exported file, the saved project,
  and closes any docked viewers.

## Scenarios

### Scenario 1: Programmatic load entry path

Steps:
1. Trigger `Bio:initBio` if not yet initialized (atlas
   `bio.lifecycle.init`). Verify `getSeqHelper()` returns the
   `ISeqHelper` singleton (atlas `bio.api.get-seq-helper`).
2. Load the FASTA file programmatically:
   `await grok.data.loadTable('System:AppData/Bio/samples/<file>.fasta')`.
   The Bio FASTA file handler (atlas `bio.io.fasta-handler`)
   ingests the file; the detector (atlas `bio.detector`) runs
   synchronously on the resulting Macromolecule column.
3. Verify detection outcome on the imported table:
   - `column.semType === 'Macromolecule'`.
   - `column.getTag('units') === 'fasta'`.
   - The FASTA cell renderer (atlas `bio.rendering.fasta`)
     paints sequence cells (no `[object Object]` fallback).

Expected:
- Programmatic load path triggers detector synchronously
  (entry-path-detector-sync invariant per atlas
  `bio.x.entry-path-detector-sync`, GROK-18616).
- Renderer dispatches based on the units tag immediately.

### Scenario 2: Drag-and-drop entry path

Steps:
1. Use Playwright's file-drop affordance to drop the same
   FASTA file onto the Datagrok window. The
   FASTA handler ingests via the drop path (atlas
   `bio.io.fasta-handler`, interaction
   `"drop .fasta file onto Datagrok"`).
2. Verify the detector classifies the Macromolecule column
   synchronously (same units / semType / renderer assertions
   as Scenario 1).
3. (Operator note) The Files-browser double-click entry path
   (a third option per atlas
   `bio.cp.fasta-import-via-multiple-entry-paths`) is
   addressable via Playwright clicking the file row in the
   Browse panel; included here as an optional verification
   when the Browse | Files surface is reachable from the test
   harness.

Expected:
- Drop entry path triggers the same synchronous detection
  outcome as the programmatic load.
- Sequence Space / MSA editor dialogs (subsequently opened)
  see a populated Macromolecule column list — the contract
  GROK-18616 reported broken on the Open-file icon path.

### Scenario 3: Export As FASTA round trip

Steps:
1. With the imported table from Scenario 1 (or 2) active,
   trigger File | Export | As FASTA... or right-click the
   Macromolecule column header → Save As FASTA... (atlas
   `bio.io.save-as-fasta`).
2. Save to the temp export path from Setup. Atlas
   `dep_lifecycle_ops[export_as_fasta]` covers this code path.
3. Re-import the exported file via
   `await grok.data.loadTable('${EXPORT_PATH}')`. The FASTA
   handler runs again on a file this scenario produced —
   round-trippable contract per atlas
   `bio.cp.import-fasta-export-fasta` and
   `bio.cp.fasta-import-via-multiple-entry-paths`.
4. Verify the re-imported table:
   - Row count matches the original FASTA's sequence count.
   - First re-imported sequence equals the first original
     sequence (string equality after `getSeqHelper`-aware
     comparison).
   - Renderer dispatches the FASTA cell renderer (atlas
     `bio.rendering.fasta`).

Expected:
- The exported FASTA is round-trippable — re-import
  reproduces the original sequences (no monomer-name drift,
  no truncation).
- atlas `dep_lifecycle_ops[import_fasta_file]` and
  `[export_as_fasta]` round-trip cleanly.

### Scenario 4: Save project with FASTA-imported table

Steps:
1. With the imported FASTA table from Scenario 1 active, save
   the project via the ribbon SAVE button (NOT Ctrl+S); project
   name from Setup; **Data Sync** toggle ON. Cancel the
   auto-share dialog if it appears. Verify
   `grok.dapi.projects.find(<id>)` returns the saved project.
   Atlas `dep_lifecycle_ops[save_project_with_analysis]`
   (`affected_source_classes: [all]`).
2. Close and reopen via
   `grok.dapi.projects.find(<id>)` → `project.open()`. Verify:
   - The Macromolecule column is back with the same
     `units` / `semType` tags (detector state survives the
     project save-reopen round trip).
   - The renderer paints FASTA cells the same as before
     (no renderer reset / no late detector mutation).

Expected:
- Project save + reopen preserves the FASTA-imported
  Macromolecule shape — the entry-path detection-sync
  invariant from Scenario 1 extends to the project-
  persistence layer.

### Scenario 5: Cleanup

Steps:
1. Delete the exported FASTA file:
   `await grok.dapi.files.delete('${EXPORT_PATH}')`.
2. Delete the project:
   `await grok.dapi.projects.delete(project)`.
3. Close any open tables / viewers / dialogs.

Expected:
- Cleanup runs in `tearDownAll` / `finally` regardless of
  earlier failures (no test-state leak across runs per the
  testing-rules in `core/docs/platform/TESTING.md`).

## Notes

- Origin: chain rev 8
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[3]`
  (`source_class: fasta_file`, `spec_target:
  bio-lifecycle-fasta-file-spec.ts`). Lands the
  intent-layer `.md` realization that
  `F-PROACTIVE-COVERAGE-01` (SCOPE_REDUCTION, cycle
  2026-06-01-bio-migrate-01) cited as the gap. `.ts`
  realization (`spec_target`) is downstream Automator /
  Phase C, not gated by F.
- atlas entry derived from
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[3]`
  (chain author, cycle 2026-05-30-bio-chain-bootstrap-01).
- `dep_lifecycle_ops` exercised (per
  `proactive_lifecycle_specs[3].bundled_ops`):
  `detect_macromolecule_on_open` (Scenario 1.3, 2.2),
  `import_fasta_file` (Scenario 1.2, 2.1, 3.3),
  `export_as_fasta` (Scenario 3.1-3.2),
  `save_project_with_analysis` (Scenario 4.1).
- target_layer rationale: `playwright` — the lifecycle spans
  drag-and-drop UI events, the Files browse path, the
  File | Export ribbon, and the Save Project dialog.
  `apitest` covers `saveAsFastaDo` round-trip
  (`fasta-export-tests.ts`) but cannot exercise the drop /
  ribbon / Browse-panel entry-path UI that GROK-18616
  surfaces.
- Sibling tests considered:
  - `public/packages/Bio/src/tests/fasta-export-tests.ts`
    (`saveAsFastaTest1`) covers `saveAsFastaDo` round-trip
    at the API layer; this scenario adds the UI entry-path
    + project-persistence layers the apitest does not
    exercise.
  - `public/packages/Bio/src/tests/detectors-tests.ts`
    covers the detector at the unit layer; this scenario
    adds the multi-entry-path detector-sync lifecycle layer.
- This scenario covers 7 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2 — well above;
  `F-STRUCT-INTERACTION-01` floor: 3 — satisfied per
  scenario). Scenario cardinality (per `## Scenarios`
  section): 5 (one cleanup + four substantive lifecycle
  scenarios) — meets the >= 2 scenarios floor.
- Manual-only subset: none of the seven covered sub_features
  appear in atlas `manual_only[]` (verified against atlas
  rev 3 `manual_only[]` list).
- `coverage_type: regression` per STEP E heuristic: this is
  general coverage of the fasta_file lifecycle shape (not a
  single critical_path golden path → not smoke; not a
  boundary value → not edge; not stress/latency-sensitive →
  not perf). Atlas critical_paths
  `bio.cp.import-fasta-export-fasta` (p1) and
  `bio.cp.fasta-import-via-multiple-entry-paths` (p1) inform
  the lifecycle steps; both have priority p1, mapping to
  `coverage_type: regression`.
- Deferrals: Files-browser double-click entry path is
  optional within Scenario 2.3 — the assertion-grade
  detector-sync verification is on the two deterministic
  entry paths (drop + programmatic). The Browse-panel
  entry path is environment-dependent for Playwright
  driving; gated on the test harness reaching that surface.
- env_requirements: `FileShare connection with per-user
  permissions` per chain `proactive_lifecycle_specs[3]`.
- Related-bug context: `related_bugs: [GROK-18616]` per
  chain `proactive_lifecycle_specs[3].bugs_reinforcing` —
  the multi-entry-path lifecycle reinforces the
  entry-path-class detection-sync gap that GROK-18616
  surfaces. Mitigated at the lifecycle layer here; the
  cross-scenario bug-focused spec
  `bio-grok-18616-spec.ts` (chain
  `bug_focused_candidates[GROK-18616]`) is the dedicated
  contract test.

---
{
  "order": 16
}
