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
    scope_reduction_proposal: |-
      One scenario step is realized via JS API substitution rather than the UI affordance the scenario step text explicitly names; the substitution has a sound technical rationale documented in the spec body header (lines 57-64) and the scenario Notes (lines 186-193), but the rationale is NOT surfaced in scenario .md frontmatter scope_reductions[] (currently []). Migrator must apply the following SR entry. SR-01 id=save-project-ribbon-jsapi-substitution check=E-TRACE-02 (S2.1) rationale=Scenario Step 2.1 names the ribbon SAVE button (NOT Ctrl+S); project name from Setup; Data Sync toggle ON; cancel the auto-share dialog if it appears. The Save Project ribbon dialog with Data Sync toggle + auto-share dialog is a platform-wide UI surface not documented in bio.md grok-browser reference selectors (verified by grep on .claude/skills/grok-browser/references/bio.md — no [name="button-Save"] / [name="dialog-Save"] / [name="dialog-Save-Project"] / [name="input-host-Data-Sync"] family); UI driving for that surface is delegated to platform-side ui-smoke scenarios elsewhere. Spec uses the canonical test-area helper helpers/projects.ts saveAllTablesWithProvenance at line 399 (mirrors sibling bio-lifecycle-macromolecule-column-spec.ts S3.3 and bio-lifecycle-fasta-file-spec.ts S4.1 per the latter's prior Critic-E SCOPE_REDUCTION at this same SR class), preserving the assertable contract (persisted project + reopened detector state per atlas dep_lifecycle_ops[save_project_with_analysis]; the [all] shorthand applies to immunum_wasm per chain proactive_lifecycle_specs[5].rationale). verdict_status=SCOPE_REDUCTION. After Migrator applies this SR entry, the spec is otherwise PASS-compliant. Mechanical checks pass: E-STRUCT-MECH-01 spec exists at expected path public/packages/UsageAnalysis/files/TestTrack/bio/bio-lifecycle-immunum-wasm-spec.ts; -02 parses cleanly (balanced braces, complete test() signature, proper imports); -03 test() block at line 103; -04 test name 'Bio immunum_wasm source-class lifecycle: init → IMGT numbering → save+reopen → re-run' substring-matches scenario heading 'Bio — immunum_wasm source-class lifecycle'; -05 imports from @playwright/test (canonical) plus the test-area-locals ../spec-login and ../helpers/projects which resolve to TestTrack/spec-login.ts and TestTrack/helpers/projects.ts — the established sibling-spec convention shared verbatim across 33+ sibling specs in TestTrack/bio/ and TestTrack/Projects/ — not an E-STRUCT-MECH-05 fabrication per the sibling Critic-E precedent (bio-lifecycle-fasta-file.md.gate-e.verdict.json + bio-lifecycle-monomer-collection.md.gate-e.verdict.json same cycle); -06 leading /* --- sub_features_covered: [...] --- */ frontmatter at lines 1-7 with 4 ids all resolving to atlas sub_features[].id entries (bio.engines.numbering-immunum bio.yaml L750, bio.annotate.numbering-scheme L640, bio.api.get-seq-helper L81, bio.lifecycle.init L71). E-TRACE-01..03: every softStep (S1.1, S1.2, S1.3-1.4, S1.4, S2.1, S2.2, S3.1-3.2, S3.3-3.4) traces to numbered scenario steps; all four Scenarios covered (Scenario 4 cleanup runs in finally per scenario Expected); verification steps map to JS API calls (expect(info.resolved).toBe(true), expect(info.result!.colCount).toBe(5), grok.dapi.projects.find post-reopen via reopenAndAssertProvenance, etc.). E-SEL-01: every [name=...] selector is class-1 in bio.md grok-browser reference per the spec-body Selector provenance block (lines 72-82) — [name="div-Bio"] L606, [name="div-Bio---Annotate"] L69, [name="div-Bio---Annotate---Apply-Numbering-Scheme..."] L421, [name="dialog-Apply-Antibody-Numbering"] L425, [name="button-OK"] standard dialog OK, [name="viewer-Grid"] standard platform selector — all verified at the cited bio.md lines. E-SEL-02: no invented selectors. E-SEL-03: no .fill() calls — all data entry surfaces go via grok.dapi.files.readCsv / grok.functions.call / dispatchEvent inside page.evaluate. E-HELP-01..02: softStep (helpers-registry L3927) and loginToDatagrok (L3933) registered; specTestOptions, stepErrors, saveAllTablesWithProvenance, reopenAndAssertProvenance, deleteProjectWithCleanup are spec-area infrastructure locals in TestTrack/spec-login.ts and TestTrack/helpers/projects.ts following the established sibling-spec convention (verified file existence at TestTrack/spec-login.ts and TestTrack/helpers/projects.ts L926, L999, L1060) — not registered-helper reinventions (no open-coded login, no open-coded project-save). E-LAYER-01..02 + E-LAYER-COMPLIANCE-01: scenario frontmatter target_layer: playwright; spec has multiple DOM-driving calls satisfying ≥1 requirement — page.locator(...).waitFor at L166, L173, L283; page.locator(...).click() at L289; page.waitForFunction at L294, L310; document.querySelector(...).click() / dispatchEvent patterns inside page.evaluate at L273-279; pyramid_layer absent in scenario frontmatter so ui-smoke sub-rule does not apply. E-BOUND-01..02: spec under allowed path public/packages/UsageAnalysis/files/TestTrack/bio/, no core/** or non-test-area writes. E-RETRY-IGNORES-GATE-B: does NOT fire — scenario frontmatter carries gate_verdicts.f (PASS) only, no gate_verdicts.b block present, therefore not a retry context per the detection predicate in spec-mode.md §'Gate B awareness on retry context'. Timeout-budget-without-reduction advisory: test.setTimeout(420_000) is above the 60-90s threshold, but the budget IS documented (spec lines 104-107) with explicit dataset-size-driven compute-wait rationale (cold Bio init ≤90s + IMGT numbering WASM compute on a tiny 47-row antibody fixture ≤30s + project save+reopen round trip + second numbering run); the dataset IS small (47 rows per spec line 113-116), so the boundary-applies note clears the soft signal — no additional SR-class advisory.
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
