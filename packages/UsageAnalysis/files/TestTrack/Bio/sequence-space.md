---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [bio.cp.sequence-space]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/sequence-space.md
migration_date: 2026-05-31
source_text_fixes:
  - bio-search-submenu-path-normalized-to-bio-analyze-per-atlas
  - sequence-space-function-label-normalized-to-sequence-space-ellipsis
  - implicit-three-dataset-matrix-made-explicit
  - dataset-filenames-changed-from-implicit-source-to-filter-prefix
  - dataset-path-pinned-to-system-appdata-bio-tests-directory
candidate_helpers: []
unresolved_ambiguities:
  - arbitrary-similarity-metric-and-method-name-edit-set-not-defined
  - sequence-space-editor-auto-column-selection-contract-not-declared
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Source Step 5 ("Change the Similarity metric and Method name
      parameters arbitrarily") is non-deterministic — atlas does not
      enumerate a canonical Similarity metric / Method name edit set for
      bio.analyze.sequence-space.editor. The edit-then-run flow is
      preserved as a step (the dialog must re-open, the inputs must be
      editable, and a second embedding result must be produced) but the
      correctness assertion on the edited-parameter result is deferred
      until atlas or operator supplies a concrete edit set; otherwise
      the verification would silently accept no-op edits.
    verdict_status: SCOPE_REDUCTION
related_bugs:
  - GROK-18616
  - GROK-19928
realized_as:
  - sequence-space-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T00:00:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T05:00:00Z
    spec_runs:
      - spec: sequence-space-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 205
        failure_keys: []
---


# Bio | Analyze | Sequence Space — sequence space integration

Covers the **Bio | Analyze | Sequence Space...** top menu across
three canonical Macromolecule notations (FASTA, HELM, MSA), with two
run modes (default parameters, then edited Similarity metric + Method
name). Verifies the multi-subsystem path: editor dialog
(`SequenceSpaceEditor`) → UMAP / t-SNE compute
(`sequenceSpaceTransform`) → ScatterPlot embedding viewer, via the
preprocess-encode engine.

The umbrella runner that exercises Sequence Space alongside Activity
Cliffs and Composition lives in `analyze.md`; the Activity-Cliffs-only
counterpart lives in `sequence-activity-cliffs.md`. This scenario is
the deeper Sequence-Space-specific scenario.

## Setup

Sequence Space is exercised on each of three test datasets to verify
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

### Scenario 1 — Open dataset and dispatch Sequence Space with defaults

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector classifies the sequence column synchronously; the table
   view opens.

2. On the menu ribbon, open **Bio** > **Analyze** > **Sequence
   Space...**. The Sequence Space dialog opens (`SequenceSpaceEditor`,
   which wraps `DimReductionBaseEditor` for the Macromolecule
   semtype). The dialog collects: Macromolecule column, Preprocessing
   function (`Encode Sequences`), Method (UMAP / t-SNE), Similarity
   metric (Hamming / Levenshtein / Monomer Chemical Distance /
   Needleman-Wunsch — or fingerprint metrics for HELM), Plot
   embeddings, Cluster embeddings.

3. Click **OK** to run with default parameters. The transform
   dispatches (`sequenceSpaceTransform`), computing embeddings via the
   preprocess-encode engine and appending embedding columns to the
   dataframe.

4. Verify the result viewer is added:
   - A ScatterPlot viewer titled "Embeddings" docks, plotting the
     2-D embedding produced by the compute.
   - Embedding columns (`Embed_X_1`, `Embed_Y_1`) are appended to
     the source dataframe.
   - When the dialog's `Cluster embeddings` was left on by default,
     a `Cluster (DBSCAN)` column is also appended.

### Scenario 2 — Re-open dialog and run with edited Similarity metric / Method name

For each dataset (same matrix as Scenario 1):

5. On the menu ribbon, re-open **Bio** > **Analyze** > **Sequence
   Space...** The dialog opens with the prior run's defaults.

6. In the dialog, change the **Similarity metric** and **Method name**
   inputs (the source text says "arbitrarily" — there's no
   canonical edit set). The inputs must remain editable (the dialog
   accepts new selections via the `<select>` widget surface) — this
   exercises the editor input re-binding contract.

7. Click **OK** to run with the edited parameters. The transform
   re-dispatches; verify a second "Embeddings" ScatterPlot docks,
   distinct from the first run's viewer, and that a second set of
   embedding columns (`Embed_X_2`, `Embed_Y_2`, and a second
   `Cluster (DBSCAN)` column if clustering is on) is appended to
   the dataframe. The edited-parameter result correctness assertion
   (e.g. that a `Levenshtein + t-SNE` run produces a meaningfully
   different embedding distribution than a `Hamming + UMAP` run) is
   deferred.

## Notes

- GROK-18616 (Sequence Space dialog showing an empty column list
  when data is loaded via the "Open file" icon) is covered more
  thoroughly by `bio-lifecycle-fasta-file.md`, which opens the same
  file through multiple entry paths.
- GROK-19928 (Sequence Space / Activity Cliffs output not surviving
  a Data-Sync project save/reopen) is covered separately in
  `bio-lifecycle-macromolecule-column.md`.
- Open questions: there's no canonical "edit set" for the Similarity
  metric / Method name re-run in Scenario 2 — the dialog is only
  checked for being editable and producing a second result, not for
  the numeric correctness of that second result. Likewise, which
  Macromolecule column the dialog auto-selects on a dataset with
  more than one Macromolecule column isn't a documented contract.

---
{
  "order": 9
}
