---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [detect_macromolecule_on_open, convert_notation]
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
---

# Bio — Macromolecule column: detect, convert, export, save & reopen

Walks a Macromolecule column through its full lifecycle: automatic
detection on open, notation conversion (FASTA → SEPARATOR), export
to FASTA and re-import, and saving/reopening a project that contains
an analysis result on the column.

Distinct from the section's analyze / convert / search scenarios
(`analyze.md`, `convert.md`, `msa.md`, etc.), which each cover a
single action end-to-end, this scenario chains several operations
together and checks that state (detected notation, converted column,
analysis output) survives save-and-reopen. It reinforces four known
bugs at the lifecycle level: renderer dispatch staying stale after a
notation convert (GROK-12164), the atomic-level/molfile conversion
path (GROK-15176), detector state going out of sync depending on how
a file was opened (GROK-18616), and analysis output not surviving a
Data-Sync project save/reopen (GROK-19928).

## Setup

- Authenticate to Datagrok as the test user (Playwright session).
- Test dataset: `System.AppData/Bio/tests/filter_FASTA.csv`
  (canonical FASTA Macromolecule column).
- Conversion target dataset: `System.AppData/Bio/tests/filter_HELM.csv`
  (HELM Macromolecule column — used as Step 2's source notation
  to exercise the HELM→SEPARATOR convert branch that hit
  GROK-12164).
- Project name: `bio-lifecycle-macromolecule-${Date.now()}`
  (unique per run to keep parallel-run safe; cleaned in Step 6).
- Service-surface dependency: `await
  grok.functions.call('Bio:getSeqHelper')` is used in Step 2 and
  Step 4 as the supported entry-point for column manipulation.
- Cleanup: Step 6 deletes the project and removes the exported
  FASTA file.

## Scenarios

### Scenario 1: Detect on open + convert-notation round trip

Steps:
1. Open `System.AppData/Bio/tests/filter_HELM.csv` via
   `grok.data.loadTable` (or Browse > Files double-click). The
   Macromolecule detector classifies the sequence column
   synchronously on open; the `units` / `aligned` / `alphabet` /
   `separator` column tags are set before any analyze dialog
   opens.
2. Verify the detector outcome:
   - `column.semType === 'Macromolecule'`.
   - `column.getTag('units') === 'helm'`.
   - The HELM cell renderer (dispatched per the detected unit)
     paints sequence cells (no `[object Object]` fallback, no
     blank canvas).
3. Right-click the Macromolecule column header > Convert ... >
   pick **SEPARATOR**. Use the dash `-` as separator — the same
   notation pair that surfaced GROK-12164.
4. Verify the convert outcome:
   - A new SEPARATOR-units Macromolecule column appears.
   - `column.getTag('units') === 'separator'`.
   - The separator cell renderer dispatches and paints correctly
     divided monomers (no monomer-collapse, no missing dividers).

Expected:
- Detector tags are stable on open (no late mutation; reading
  `units` immediately after load returns the final value).
- Convert mutates the column shape and the renderer-dispatch
  follows the new tags consistently — the
  detector-renderer-after-convert contract (GROK-12164) is
  satisfied at the lifecycle layer.

### Scenario 2: Import FASTA → Export FASTA → re-import round trip

Steps:
1. Open `System.AppData/Bio/tests/filter_FASTA.csv` (canonical
   FASTA Macromolecule sample). The Bio FASTA file handler
   ingests; the detector classifies synchronously.
2. Right-click the Macromolecule column > Save As FASTA... (or
   File > Export > As FASTA...). Save to a temp path (e.g.
   `System.AppData/UsageAnalysis/temp/lifecycle-${Date.now()}.fasta`).
3. Re-open the exported FASTA file via `grok.data.loadTable`
   against the temp path. The FASTA handler runs again on a file
   that this scenario itself produced — a round-trippable
   contract.
4. Verify the imported table:
   - Row count matches the original FASTA's sequence count
     (modulo any silent dedupe — assert exact equality).
   - The first re-imported sequence, normalized, equals the
     first original sequence (string equality after
     `getSeqHelper`-aware comparison).
   - The Macromolecule renderer dispatches the FASTA cell
     renderer.

Expected:
- The exported FASTA is round-trippable: re-importing reproduces
  the original sequence content (no monomer-name drift, no
  separator artifacts, no truncation).
- The entry-path detector-sync invariant holds for the
  programmatic load path (`grok.data.loadTable`) (GROK-18616) —
  the dialog side is exercised in Scenario 3 via Sequence Space.

### Scenario 3: Save project with analysis + reopen restores analysis output

Steps:
1. With the FASTA table from Scenario 2 active (or re-open
   `filter_FASTA.csv` if Scenario 2 was not run), trigger
   **Bio | Analyze | Sequence Space...**. Run with default
   parameters. Wait for the embedding compute to finish; the
   ScatterPlot embedding viewer docks.
2. Verify the embedding output is present:
   - The two new embedding columns (`Embed_X`, `Embed_Y` or
     equivalent) are added to the DataFrame.
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
  analysis output shape (GROK-19928) at the lifecycle layer.
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

- Sibling coverage: `fasta-export-tests.ts` covers the FASTA export
  round trip at the API layer, `sequence-space-tests.ts` covers the
  embedding-compute API, and `detectors-tests.ts` / `convert-tests.ts`
  cover the detector and notation-convert logic at the unit layer.
  This scenario adds the UI-driven, cross-op chaining, and
  save-and-reopen persistence layers those tests don't exercise.
- This scenario doesn't cover the macromolecule-difference (diff-cell)
  renderer, even though the convert-notation invariant is adjacent to
  GROK-16596 reset-state behavior — that renderer's correctness is a
  manual/visual check only.
- Deferrals: reopening the project is done via the JS API
  (`project.open()`) rather than clicking through the Browse >
  Dashboards UI — the full UI-driven reopen flow is covered by
  ui-smoke scenarios elsewhere in the section.
- The full cross-package PubChem-standardization contract for
  GROK-15176 (isotope-flag validity on generated molfiles) is not
  covered by a dedicated spec; this scenario only reinforces the
  invariant at the convert / round-trip persistence layer.

---
{
  "order": 13
}
