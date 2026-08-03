---
feature: bio
sub_features_covered:
  - bio.engines.numbering-immunum
  - bio.annotate.numbering-scheme
  - bio.api.get-seq-helper
  - bio.lifecycle.init
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-lifecycle-immunum-wasm-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T14:30:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T06:05:00Z
    spec_runs:
      - spec: bio-lifecycle-immunum-wasm-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 144
        failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T13:15:00Z
    failure_keys: []
---

# Bio — `immunum_wasm` source-class lifecycle

Proactive lifecycle scenario for the `immunum_wasm` source
class per
`scenario-chains/bio.yaml#proactive_lifecycle_specs[5]`. Walks
the immunum WASM module through the only non-agnostic
`dep_lifecycle_op` declared in atlas
`dep_lifecycle_ops[].affected_source_classes` that touches the
shape (the `[all]` shorthand on
`save_project_with_analysis`): load WASM in worker (bundled
with package version) → IMGT/Kabat numbering on an antibody
column via the Immunum engine
(`bio.engines.numbering-immunum`) → save project with the
numbering output → reopen and verify the numbering annotation
survives the round trip.

`external_deps: []` per chain `proactive_lifecycle_specs[5]`
— the WASM asset travels with the package version and runs
in-process in a dedicated worker; no runtime entity-type
dependency. `env_requirements: []` — runs against a synthetic
antibody-sequences dataset.

This is the smallest of the six proactive lifecycle cells per
chain rationale: "Only one non-agnostic op affects this source
class (`save_project_with_analysis` via [all] shorthand) since
the WASM asset travels with the package version and has no
runtime entity-type dep — load is in-process."

## Setup

- Authenticate to Datagrok as the test user (Playwright session).
- Test dataset: a small antibody-sequences DataFrame. Inline
  construct via `DG.DataFrame.fromCsv(...)` with a tiny
  hand-built antibody-sequences CSV (one or two heavy-chain
  rows, FASTA notation) if no AppData antibody sample is
  available. The Immunum engine accepts any valid antibody
  sequence input per atlas
  `bio.engines.numbering-immunum`.
- Project name:
  `bio-lifecycle-immunum-wasm-${Date.now()}`.
- Service-surface dependency:
  `await grok.functions.call('Bio:getSeqHelper')` is used in
  Step 1 (atlas `bio.api.get-seq-helper`).
- Cleanup: Step 4 deletes the saved project and closes any
  open viewers.

## Scenarios

### Scenario 1: Load WASM and run numbering

Steps:
1. Trigger `Bio:initBio` if not yet initialized (atlas
   `bio.lifecycle.init`). Verify `getSeqHelper()` resolves
   (atlas `bio.api.get-seq-helper`).
2. Materialize the antibody-sequences DataFrame in memory
   (CSV-from-string into `DG.DataFrame`). Verify the
   detector classifies the sequence column as Macromolecule
   (atlas `bio.detector` — not in this scenario's
   `sub_features_covered` because the immunum lifecycle is
   the primary surface, not the detector).
3. Trigger ribbon `Bio | Annotate | Apply Numbering
   Scheme...` (atlas `bio.annotate.numbering-scheme`). In the
   dialog, select the Immunum engine (`meta.role:
   antibodyNumbering` registration; per atlas
   `bio.engines.numbering-immunum`) and the IMGT scheme.
4. Verify the numbering result. Per atlas
   `bio.engines.numbering-immunum`:
   - A 5-column result DataFrame is produced
     (`position_names`, `chain_type`, `annotations_json`,
     `numbering_detail`, `numbering_map`).
   - All five columns are populated with non-null values for
     each input antibody row.
   - The WASM worker call completes within a reasonable
     time bound (immunum runs in-process per the worker
     model — no Docker container dependency).

Expected:
- The Immunum WASM module loads on first call (lazy init
  per the WASM-loader pattern in
  `public/.claude/rules/wasm.md`).
- IMGT numbering output shape matches the atlas-declared
  5-column contract.

### Scenario 2: Save project with numbering output

Steps:
1. With the antibody table + numbering output from Scenario
   1 active, save the project via the ribbon SAVE button
   (NOT Ctrl+S); project name from Setup; **Data Sync**
   toggle ON. Cancel the auto-share dialog if it appears.
   Verify `grok.dapi.projects.find(<id>)` returns the saved
   project. Atlas
   `dep_lifecycle_ops[save_project_with_analysis]`
   (`affected_source_classes: [all]` shorthand applies to
   `immunum_wasm` per chain
   `proactive_lifecycle_specs[5].rationale`).
2. Verify the project record persists the numbering output
   columns (column shape + values).

Expected:
- Project save succeeds with the Immunum numbering output
  intact.

### Scenario 3: Reopen project + WASM re-load

Steps:
1. Close the project. Reopen via
   `grok.dapi.projects.find(<id>)` → `project.open()`.
2. Verify the reopened project state:
   - The antibody-sequences table is back.
   - The 5-column numbering output is back with the same
     values as before save.
   - The Macromolecule column tags are intact.
3. Run the Immunum numbering again on the reopened table to
   verify the WASM module can re-load cleanly in the new
   session (the bundled-with-package-version WASM asset
   loads via the same code path on every fresh init).
4. Verify the re-run produces the same numbering output as
   the original Scenario 1 run (deterministic).

Expected:
- Project reopen restores the numbering output; the WASM
  re-load on a fresh session is deterministic and produces
  identical numbering for identical input.

### Scenario 4: Cleanup

Steps:
1. Delete the project:
   `await grok.dapi.projects.delete(project)`.
2. Close any open tables / viewers.

Expected:
- Cleanup runs in `tearDownAll` / `finally` regardless of
  earlier failures (no test-state leak across runs per the
  testing-rules in `core/docs/platform/TESTING.md`).

## Notes

- Origin: chain rev 8
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[5]`
  (`source_class: immunum_wasm`, `spec_target:
  bio-lifecycle-immunum-wasm-spec.ts`). Lands the
  intent-layer `.md` realization that
  `F-PROACTIVE-COVERAGE-01` (SCOPE_REDUCTION, cycle
  2026-06-01-bio-migrate-01) cited as the gap. `.ts`
  realization (`spec_target`) is downstream Automator /
  Phase C, not gated by F.
- atlas entry derived from
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[5]`
  (chain author, cycle 2026-05-30-bio-chain-bootstrap-01).
- `dep_lifecycle_ops` exercised (per
  `proactive_lifecycle_specs[5].bundled_ops`):
  `save_project_with_analysis` (Scenario 2.1). Only one
  non-agnostic op affects this source class — chain
  rationale clarifies the lifecycle surface is minimal
  because the WASM asset is bundled with the package
  version and load is in-process (no separate
  load_xxx / save_xxx op against immunum_wasm).
- target_layer rationale: `playwright` — the lifecycle
  spans the `Bio | Annotate | Apply Numbering Scheme...`
  ribbon path + the Save Project ribbon path; both
  surfaces require DOM driving. `apitest` covers the
  Immunum engine result-shape at the API layer but does
  not exercise the ribbon dispatch nor the project save
  UI. JS API substitutes are used for project persistence
  assertions per the same pattern as sibling
  `bio-lifecycle-*.md` scenarios.
- Sibling tests considered:
  - Antibody numbering engine API tests (if present under
    `public/packages/Bio/src/tests/`) cover the
    `applyNumberingScheme` engine surface at the API
    layer; this scenario adds the UI ribbon-dispatch +
    project-persistence lifecycle layer.
- This scenario covers 4 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2 — above;
  `F-STRUCT-INTERACTION-01` floor: 3 — satisfied per
  scenario). Scenario cardinality (per `## Scenarios`
  section): 4 (one cleanup + three substantive lifecycle
  scenarios) — meets the >= 2 scenarios floor.
- Manual-only subset: none of the four covered sub_features
  appear in atlas `manual_only[]` (verified against atlas
  rev 3 `manual_only[]` list).
- `coverage_type: regression` per STEP E heuristic: this is
  general coverage of the immunum_wasm lifecycle shape (not
  a single critical_path golden path → not smoke; not a
  boundary value → not edge; not stress/latency-sensitive
  → not perf). Atlas `bio.cp.numbering-scheme` (priority
  p1) informs the workflow; priority p1 maps to
  `coverage_type: regression`.
- Deferrals: none mandatory. The WASM-loading mechanism
  (lazy init pattern) is exercised implicitly via the
  first-call path in Scenario 1.3 and the re-load path in
  Scenario 3.3 — no explicit WASM module-handle assertions
  are made beyond observable outputs (would require
  test-internal hooks).
- env_requirements: `[]` per chain
  `proactive_lifecycle_specs[5]`. No external dependencies
  — WASM asset travels with the package version, loads
  in-process in a dedicated worker.
- Related-bug context: `related_bugs: []` per chain
  `proactive_lifecycle_specs[5].bugs_reinforcing` (empty).

---
{
  "order": 18
}
