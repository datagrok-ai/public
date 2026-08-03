---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [bio.cp.activity-cliffs]
realizes: [bio.analyze.activity-cliffs]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/sequence-activity-cliffs.md
migration_date: 2026-05-31
source_text_fixes:
  - bio-search-submenu-path-normalized-to-bio-analyze-per-atlas
  - sequence-activity-cliffs-function-label-normalized-to-activity-cliffs-ellipsis
  - implicit-three-dataset-matrix-made-explicit
  - dataset-filenames-changed-from-sample-prefix-to-filter-prefix
  - dataset-path-pinned-to-system-appdata-bio-tests-directory
candidate_helpers: []
unresolved_ambiguities:
  - arbitrary-similarity-and-method-name-edit-set-not-defined
  - activities-column-default-may-pick-length-over-activity
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Source Step 5 ("Change the Similarity and Method name parameters
      arbitrarily") is non-deterministic — atlas does not enumerate a
      canonical Similarity / Method edit set for
      bio.analyze.activity-cliffs.editor. The edit-then-run flow is
      preserved as a step (the dialog must re-open, the inputs must be
      editable, and a second cliff result must be produced) but the
      correctness assertion on the edited-parameter result is deferred
      until atlas or operator supplies a concrete edit set; otherwise the
      verification would silently accept no-op edits.
    verdict_status: SCOPE_REDUCTION
related_bugs:
  - GROK-19150
  - GROK-19928
  - GROK-16111
realized_as:
  - sequence-activity-cliffs-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T13:30:00Z
    failure_keys: []
  a:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T14:00:00Z
    failure_keys: []
    review_round: 1
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T07:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T08:46:00Z
    spec_runs:
      - spec: sequence-activity-cliffs-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 231
        failure_keys: []
---

# Bio | Analyze | Activity Cliffs — sequence activity cliffs integration

Covers the **Bio | Analyze | Activity Cliffs...** top menu across
three canonical Macromolecule notations (FASTA, HELM, MSA), with two
run modes (default parameters, then edited Similarity + Method name).
Verifies the multi-subsystem path: editor dialog → embedding compute
→ ScatterPlot viewer with cliff overlay.

The umbrella runner that exercises Activity Cliffs alongside Sequence
Space and Composition lives in `analyze.md`; the Sequence-Space-only
counterpart lives in `sequence-space.md`. This scenario is the
deeper Activity-Cliffs-specific scenario.

## Setup

Activity Cliffs is exercised on each of three test datasets to verify
notation-agnostic behavior across the canonical Macromolecule
notations:

- `System.AppData/Bio/tests/filter_FASTA.csv`
- `System.AppData/Bio/tests/filter_HELM.csv`
- `System.AppData/Bio/tests/filter_MSA.csv`

The implicit matrix is **3 datasets × 2 run modes (defaults, edited
parameters) = 6 cells**. The defaults run is verified per cell; the
edited-parameters run is preserved as a flow but its correctness
assertion is deferred.

## Scenarios

### Scenario 1 — Open dataset and dispatch Activity Cliffs with defaults

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector classifies the sequence column synchronously; the table
   view opens.

2. On the menu ribbon, open **Bio** > **Analyze** > **Activity
   Cliffs...**. The Activity Cliffs dialog opens
   (`SeqActivityCliffsEditor`). The dialog collects: Macromolecule
   column, Activities (numeric), Similarity metric, Method,
   similarity cutoff.

3. Click **OK** to run with default parameters. The transform
   dispatches (`seqActivityCliffsTransform`), computing embeddings and
   stamping the `seqActivityCliffsParams` tag on the table.

4. Verify the result viewer is added:
   - A ScatterPlot viewer titled "Activity cliffs" docks
     (`seqActivityCliffsInitFunction`) with the cliff overlay drawn
     on top of the sequence-space embedding.
   - Embedding columns (X / Y) and the SALI scoring column are
     appended to the source dataframe.

### Scenario 2 — Re-open dialog and run with edited Similarity / Method

For each dataset (same matrix as Scenario 1):

5. On the menu ribbon, re-open **Bio** > **Analyze** > **Activity
   Cliffs...** The dialog opens with the prior run's defaults.

6. In the dialog, change the **Similarity** metric and **Method name**
   inputs (the source text says "arbitrarily" — there's no canonical
   edit set). The inputs must remain editable (the dialog accepts new
   selections via the `<select>` widget surface) — this exercises the
   editor input re-binding contract.

7. Click **OK** to run with the edited parameters. The transform
   re-dispatches; verify a second "Activity cliffs" ScatterPlot
   docks, distinct from the first run's viewer (total Activity Cliffs
   viewer count = 2). The edited-parameter result correctness
   assertion (e.g. that a `Levenshtein + t-SNE` run produces a
   meaningfully different SALI distribution than a `Hamming + UMAP`
   run) is deferred.

## Notes

- GROK-19150 (Activity Cliffs multi-instance click-routing — clicking
  on a second analysis routes to the first analysis's results grid)
  is not covered here; this scenario only opens one Activity Cliffs
  analysis at a time. A dedicated regression spec for it is planned.
- GROK-19928 (Sequence Space / Activity Cliffs output not surviving a
  Data-Sync project save/reopen) is covered separately in
  `bio-lifecycle-macromolecule-column.md`.
- GROK-16111 (empty/null current-row input not rejected with an
  error) is covered separately in `empty-input-row-viewers.md`; this
  scenario only uses non-empty fixtures.
- Open questions: there's no canonical "edit set" for the Similarity
  metric / Method name re-run in Scenario 2 — the dialog is only
  checked for being editable and producing a second result, not for
  the numeric correctness of that second result. Likewise, which
  column the dialog's Activities input defaults to on a
  multi-numeric-column dataset isn't a documented contract (on the
  current test fixtures it happens to pick `Length` rather than
  `Activity`).

---
{
  "order": 8
}
