---
feature: bio
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [bio.cp.sequence-space, bio.cp.activity-cliffs, bio.cp.composition-analysis]
realizes: [bio.analyze.sequence-space, bio.analyze.activity-cliffs, bio.analyze.composition, bio.weblogo]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/analyze.md
migration_date: 2026-05-31
source_text_fixes:
  - duplicated-step-2-renumbered-to-3-and-4
  - implicit-three-dataset-matrix-made-explicit
  - bio-analyze-submenu-path-normalized-to-atlas-form
candidate_helpers: []
realized_as:
  - analyze-spec.ts
unresolved_ambiguities:
  - arbitrary-changed-parameters-edit-set-not-defined
  - composition-context-panel-settings-checklist-not-defined
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Source Step 3 ("Run the function once more with arbitrary changed
      parameters") is non-deterministic — atlas does not enumerate a
      canonical parameter-edit set for Sequence Space, Activity Cliffs,
      or Composition. Deferred until atlas or operator supplies a
      concrete edit set; otherwise the verification would silently
      accept no-op edits.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: A-CONT-01
    rationale: |
      Source Step 2 (second instance, renumbered to Step 4 in migrated
      form) — "check all settings on the Context Panel" — is open-ended;
      atlas does not enumerate which Context-Pane property surface the
      Composition viewer Gear icon exposes. Deferred until atlas or
      operator supplies a concrete checklist of Context-Pane properties
      to verify.
    verdict_status: SCOPE_REDUCTION
related_bugs:
  - GROK-18616
  - GROK-19928
  - GROK-19150
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T12:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T13:00:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T12:11:17Z
    spec_runs:
      - spec: analyze-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 123
        failure_keys: []
---

# Bio | Analyze — umbrella integration smoke

Umbrella integration scenario for the three `Bio | Analyze` top-menu
functions (Sequence Space, Activity Cliffs, Composition) across three
canonical Macromolecule notations (FASTA, HELM, MSA). Verifies
multi-subsystem co-existence: top-menu dispatch → editor dialogs →
viewer docking → Context-Pane wiring (Composition Gear).

The deeper per-function scenarios live in `sequence-space.md`,
`sequence-activity-cliffs.md`, and `composition-analysis.md`; this
scenario is the cross-function runner.

## Setup

All functions are exercised on each of:

- `System.AppData/Bio/tests/filter_FASTA.csv`
- `System.AppData/Bio/tests/filter_HELM.csv`
- `System.AppData/Bio/tests/filter_MSA.csv`

The implicit matrix is **3 datasets × 3 analyze functions = 9 cells**.
Per cell, the run-with-defaults assertion is verified; the
run-with-edited-parameters assertion is deferred.

## Scenarios

### Scenario 1 — Open dataset and dispatch each Analyze function

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector classifies the sequence column synchronously (atlas
   `bio.detector`); the table view opens.

2. On the menu ribbon, open **Bio** > **Analyze**. Verify the submenu
   exposes:
   - **Sequence Space...**
   - **Activity Cliffs...**
   - **Composition**

### Scenario 2 — Run each function with default parameters

For each dataset × each of the three Analyze functions:

3. Invoke the function from **Bio** > **Analyze**. For **Sequence
   Space...** and **Activity Cliffs...**, a dialog opens; click **OK**
   to run with default parameters. **Composition** docks its viewer
   without a prompt-style dialog.

4. Verify a viewer opens / docks:
   - Sequence Space → ScatterPlot with embedding columns appended.
   - Activity Cliffs → ScatterPlot with cliff overlay.
   - Composition → WebLogo viewer docked (viewer presence is
     verified; pixel-level WebLogo paint is a manual/visual check and
     not asserted here).

### Scenario 3 — Composition Gear → Context Panel wiring

5. After running **Composition** in Scenario 2 on `filter_FASTA.csv`,
   on the docked Composition viewer click the **Gear** icon. Verify the
   Context Panel opens with the Composition viewer's property surface
   bound (presence of the property editor surface only; the concrete
   property checklist is deferred).

## Notes

- **Manual-only surface touched but not asserted in pixel form:** the
  WebLogo viewer painted by Composition. Its letter/color/conservation
  paint is human-inspection only; this scenario asserts viewer-docking
  presence and the Gear → Context Panel wiring only, not pixel-level
  correctness.
- **Deliberately not verified:** re-running each function a second
  time with arbitrary changed parameters (there's no fixed set of
  "changed" values to check results against), and the exact list of
  properties shown on the Composition viewer's Context Panel after
  clicking the Gear icon. Both flows still run end-to-end (the dialog
  reopens, the Context Panel opens); only the fine-grained correctness
  assertion on the edited state is deferred, pending a concrete
  checklist to verify against.

---
{
  "order": 1
}
