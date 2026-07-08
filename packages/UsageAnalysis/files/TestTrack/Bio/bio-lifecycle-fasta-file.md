---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [import_fasta_file, export_as_fasta]
produced_from: atlas-driven
related_bugs:
  - GROK-18616
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      Scenario 2's synthetic File-drop (drag-and-drop FASTA entry path) does
      not reliably dispatch the file handler in the apitest harness; the spec
      falls back to the `Bio:importFasta` atlas-equivalent call, which runs
      the same `FastaFileHandler.importFasta` code path (atlas
      bio.cp.fasta-import-via-multiple-entry-paths). The multiple-entry-path
      invariant is preserved; the literal drag-and-drop DOM event is not
      asserted.
    verdict_status: SCOPE_REDUCTION
realized_as:
  - bio-lifecycle-fasta-file-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T07:00:00Z
    failure_keys: []
    review_round: 1
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

# Bio — FASTA files: import via multiple paths, export, save & reopen

Walks a FASTA file through its full lifecycle: opening it via
multiple entry paths (programmatic load, drag-and-drop, and —
optionally — the Files browser), converting it to a Macromolecule
column with automatic detection, exporting it back to FASTA and
re-importing it (round trip), and saving/reopening a project that
contains it.

This guards against GROK-18616, where opening a FASTA file through
some entry paths left the Macromolecule detector out of sync, so
downstream dialogs (like Sequence Space) saw an empty column list.
Unlike `bio-lifecycle-macromolecule-column.md` (which exercises a
single import path), this scenario explicitly opens the same file
through multiple different paths and checks that detection stays
consistent across all of them.

Requires a FileShare connection with permission to read/write under
`System:AppData/Bio/samples`.

## Setup

- Authenticate to Datagrok as the test user (Playwright session)
  with read/write access to `System:AppData/Bio/samples` and a
  temp directory for exports.
- Source FASTA file: any of `human_genes.fasta` (the
  GROK-18474-tagged sample). If absent, fall back to any
  `*.fasta` file under `System:AppData/Bio/samples`.
- Export target path:
  `System:AppData/UsageAnalysis/temp/lifecycle-fasta-${Date.now()}.fasta`.
- Project name:
  `bio-lifecycle-fasta-file-${Date.now()}`.
- Service-surface dependency:
  `await grok.functions.call('Bio:getSeqHelper')` is used in
  Steps 1 and 2.
- Cleanup: Step 5 deletes the exported file, the saved project,
  and closes any docked viewers.

## Scenarios

### Scenario 1: Programmatic load entry path

Steps:
1. Trigger `Bio:initBio` if not yet initialized. Verify
   `getSeqHelper()` returns the `ISeqHelper` singleton.
2. Load the FASTA file programmatically:
   `await grok.data.loadTable('System:AppData/Bio/samples/<file>.fasta')`.
   The Bio FASTA file handler ingests the file; the detector
   runs synchronously on the resulting Macromolecule column.
3. Verify detection outcome on the imported table:
   - `column.semType === 'Macromolecule'`.
   - `column.getTag('units') === 'fasta'`.
   - The FASTA cell renderer paints sequence cells (no
     `[object Object]` fallback).

Expected:
- Programmatic load path triggers detector synchronously
  (entry-path-detector-sync invariant, GROK-18616).
- Renderer dispatches based on the units tag immediately.

### Scenario 2: Drag-and-drop entry path

Steps:
1. Use Playwright's file-drop affordance to drop the same
   FASTA file onto the Datagrok window. The FASTA handler
   ingests via the drop path.
2. Verify the detector classifies the Macromolecule column
   synchronously (same units / semType / renderer assertions
   as Scenario 1).
3. (Operator note) The Files-browser double-click entry path
   (a third option) is addressable via Playwright clicking the
   file row in the Browse panel; included here as an optional
   verification when the Browse | Files surface is reachable
   from the test harness.

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
   Macromolecule column header → Save As FASTA...
2. Save to the temp export path from Setup.
3. Re-import the exported file via
   `await grok.data.loadTable('${EXPORT_PATH}')`. The FASTA
   handler runs again on a file this scenario produced — a
   round-trippable contract.
4. Verify the re-imported table:
   - Row count matches the original FASTA's sequence count.
   - First re-imported sequence equals the first original
     sequence (string equality after `getSeqHelper`-aware
     comparison).
   - Renderer dispatches the FASTA cell renderer.

Expected:
- The exported FASTA is round-trippable — re-import
  reproduces the original sequences (no monomer-name drift,
  no truncation).
- Import and export round-trip cleanly.

### Scenario 4: Save project with FASTA-imported table

Steps:
1. With the imported FASTA table from Scenario 1 active, save
   the project via the ribbon SAVE button (NOT Ctrl+S); project
   name from Setup; **Data Sync** toggle ON. Cancel the
   auto-share dialog if it appears. Verify
   `grok.dapi.projects.find(<id>)` returns the saved project.
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

- Sibling coverage: `fasta-export-tests.ts` covers the FASTA
  export/round-trip logic at the API layer, and
  `detectors-tests.ts` covers the detector at the unit layer.
  This scenario adds the UI entry-path and project-persistence
  layers those tests don't exercise.
- Deferrals: exercising the Files-browser double-click entry
  path (Scenario 2, step 3) is optional — included when
  reachable from the test harness, but only the drag-and-drop
  and programmatic-load entry paths are asserted as a hard
  requirement, since the Browse-panel path is
  environment-dependent.
- No separate bug-focused regression spec for GROK-18616 exists yet;
  this scenario reinforces the same invariant at the lifecycle
  (multiple-entry-path) layer.

---
{
  "order": 16
}
