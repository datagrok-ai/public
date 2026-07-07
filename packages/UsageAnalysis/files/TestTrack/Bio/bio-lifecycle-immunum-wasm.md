---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [save_project_with_analysis]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      This scenario does not assert on the WASM module's internal state
      (e.g. whether a fresh session re-uses vs. re-loads the numbering-engine
      module handle) — only on its observable outputs (the numbering result).
      Asserting internal module state would require test-only hooks, so it is
      deferred.
    verdict_status: SCOPE_REDUCTION
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

# Bio — Antibody numbering (Immunum WASM): run, save & reopen

Runs antibody numbering (IMGT/Kabat schemes) via the Immunum WASM
engine on an antibody sequence column, then checks that the
numbering output survives a project save and reopen.

The Immunum WASM module ships bundled with the Bio package and loads
in-process in a worker — no external server or file-share dependency
is needed, so this scenario runs against a small synthetic antibody
dataset built in-memory.

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
  Step 1.
- Cleanup: Step 4 deletes the saved project and closes any
  open viewers.

## Scenarios

### Scenario 1: Load WASM and run numbering

Steps:
1. Trigger `Bio:initBio` if not yet initialized. Verify
   `getSeqHelper()` resolves.
2. Materialize the antibody-sequences DataFrame in memory
   (CSV-from-string into `DG.DataFrame`). Verify the
   detector classifies the sequence column as Macromolecule.
3. Trigger ribbon `Bio | Annotate | Apply Numbering
   Scheme...`. In the dialog, select the Immunum engine
   (`meta.role: antibodyNumbering` registration) and the IMGT
   scheme.
4. Verify the numbering result:
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
   project.
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

- Sibling coverage: dedicated API-level tests for the Immunum
  numbering engine (if present under `Bio/src/tests/`) cover the
  `applyNumberingScheme` call directly; this scenario adds the UI
  ribbon-dispatch and project-persistence layers those tests don't
  exercise.
- Deferrals: this scenario doesn't assert on the WASM module's
  internal state (e.g. whether a fresh session re-uses vs. re-loads
  the module handle) — only on its observable outputs (the numbering
  result). Asserting internal module state would require test-only
  hooks.

---
{
  "order": 18
}
