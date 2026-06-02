---
feature: bio
sub_features_covered:
  - bio.analyze.sequence-space
  - bio.analyze.sequence-space.top-menu
  - bio.analyze.sequence-space.editor
  - bio.analyze.sequence-space.transform
target_layer: playwright
coverage_type: regression
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

Integration scenario for the `Bio | Analyze | Sequence Space...` top
menu (atlas `bio.analyze.sequence-space.top-menu` —
`package.ts#L740`) across three canonical Macromolecule notations
(FASTA, HELM, MSA), with two run modes (default parameters, then
edited Similarity metric + Method name). Verifies the multi-subsystem
path: editor dialog (`bio.analyze.sequence-space.editor`,
`SequenceSpaceEditor` — `package.ts#L237`) → UMAP / t-SNE compute
(`bio.analyze.sequence-space.transform`, `sequenceSpaceTransform` —
`package.ts#L789`) → ScatterPlot embedding viewer → preprocess-encode
engine (`bio.engines.preprocess-encode`).

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
assertion is deferred (see SR-01).

## Scenarios

### Scenario 1 — Open dataset and dispatch Sequence Space with defaults

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector (atlas `bio.detector`) classifies the sequence column
   synchronously; the table view opens.

2. On the menu ribbon, open **Bio** > **Analyze** > **Sequence
   Space...**. The Sequence Space dialog opens (atlas
   `bio.analyze.sequence-space.editor`, `SequenceSpaceEditor` —
   `package.ts#L237`, which wraps `DimReductionBaseEditor` for the
   Macromolecule semtype). The dialog collects: Macromolecule column,
   Preprocessing function (`Encode Sequences`), Method (UMAP / t-SNE),
   Similarity metric (Hamming / Levenshtein / Monomer Chemical
   Distance / Needleman-Wunsch — or fingerprint metrics for HELM),
   Plot embeddings, Cluster embeddings.

3. Click **OK** to run with default parameters. The transform
   dispatches (atlas `bio.analyze.sequence-space.transform`,
   `sequenceSpaceTransform` — `package.ts#L789`), computing
   embeddings via the preprocess-encode engine
   (`bio.engines.preprocess-encode`) and appending embedding columns
   to the dataframe.

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
   inputs (the source text says "arbitrarily" — atlas does not pin
   a canonical edit set, see SR-01 and `unresolved_ambiguities`).
   The inputs must remain editable (the dialog accepts new
   selections via the `<select>` widget surface) — this exercises
   the editor input re-binding contract.

7. Click **OK** to run with the edited parameters. The transform
   re-dispatches; verify a second "Embeddings" ScatterPlot docks,
   distinct from the first run's viewer, and that a second set of
   embedding columns (`Embed_X_2`, `Embed_Y_2`, and a second
   `Cluster (DBSCAN)` column if clustering is on) is appended to
   the dataframe. The edited-parameter result correctness assertion
   (e.g. that a `Levenshtein + t-SNE` run produces a meaningfully
   different embedding distribution than a `Hamming + UMAP` run) is
   deferred — see SR-01.

## Notes

- **Source-text fixes silently applied** during migration (recorded
  in frontmatter `source_text_fixes`):
  - The original scenario placed Sequence Space under **Bio > Search
    > Sequence Space**. Atlas registers this surface as
    `bio.analyze.sequence-space.top-menu` under **Bio | Analyze |
    Sequence Space...** (`package.ts#L740`); the prior run log
    (`sequence-space-run.md`) confirms that the menu path "Bio >
    Search" does not exist on the live platform and that the live
    function lives at "Bio > Analyze > Sequence Space...". The
    chain flagged this as an unresolved ambiguity for operator
    confirmation; atlas + run-log evidence converge, so the migrated
    scenario uses the atlas-canonical path.
  - The original scenario referred to the function as "Sequence
    Space" (no ellipsis). The atlas-canonical menu label is
    "Sequence Space..." with the trailing ellipsis denoting a
    dialog-opening function. The migrated scenario uses
    "Sequence Space..." consistent with atlas
    `bio.analyze.sequence-space.top-menu` (`package.ts#L740`).
  - The implicit 3-dataset matrix was listed in Step 1 ("Test for
    data: ...") but not explicitly framed as a per-dataset loop;
    surfaced here under Setup as a 3 × 2 matrix.
  - The original scenario named the datasets `filter_FASTA.csv`,
    `filter_HELM.csv`, `filter_MSA.csv` without a directory path.
    The migrated scenario pins them to `System.AppData/Bio/tests/`
    to match the `filter_*.csv` family resolution pattern used by
    the sister `analyze-spec.ts`, `sequence-space-spec.ts`, and
    `sequence-activity-cliffs-spec.ts`. The prior run log
    (`sequence-space-run.md`) used the `Bio/samples/` family
    (`FASTA.csv` / `HELM.csv` / `MSA.csv`) instead — a different
    dataset family; the divergence between run-log empirical
    resolution and migrated fixture choice is acknowledged here.
    Operator clarification is welcomed if the `samples/` family is
    preferred over the purpose-built `tests/` family for the
    migrated scenario.
  - The original scenario gave only bare filenames with no
    directory path. The migrated scenario pins the path to
    `System.AppData/Bio/tests/` to match the `filter_*.csv` family
    selection above and to align with the sister specs' resolution
    pattern.
- **Sub-features covered:**
  - `bio.analyze.sequence-space` (atlas L388) — feature root.
  - `bio.analyze.sequence-space.top-menu` (atlas L397) — `Bio |
    Analyze | Sequence Space...` dispatch (`package.ts#L740`).
  - `bio.analyze.sequence-space.editor` (atlas L414) —
    `SequenceSpaceEditor` dialog wrapping `DimReductionBaseEditor`
    (`package.ts#L237`).
  - `bio.analyze.sequence-space.transform` (atlas L406) —
    `sequenceSpaceTransform` embedding compute (`package.ts#L789`).
  Maps directly onto atlas critical path `bio.cp.sequence-space`
  (p0, `derived_from: public/packages/Bio/src/package.ts#L740`),
  which also references `bio.engines.preprocess-encode` (exercised
  transitively when the transform runs).
- **Related bugs** surfaced for downstream awareness; per-bug
  cross-cutting invariants are delegated to chain-level
  `bug_focused_candidates`:
  - GROK-18616 — Sequence Space dialog opens with an empty column
    list when data is loaded via the "Open file" icon (entry-path
    detector-sync gap). This scenario opens datasets through a
    single entry path (the `System.AppData/Bio/tests/` Files-style
    load); the multi-entry-path invariant requires opening the same
    dataset via Files browser, Open-file icon, drag-and-drop, and
    programmatic `grok.data.loadTable` and is delegated to
    `bio-grok-18616-spec.ts` (chain
    `bug_focused_candidates[GROK-18616]` spans
    `analyze.md:Step 1` + `sequence-space.md:Step 2`). Atlas
    cross-cutting: `bio.x.entry-path-detector-sync`; atlas critical
    path: `bio.cp.fasta-import-via-multiple-entry-paths` (p1,
    `derived_from: bug-library/bio.yaml#GROK-18616`).
  - GROK-19928 — Sequence Space and Activity Cliffs analysis output
    not properly saved/restored in datasync projects (reopen crashes
    with `NullError` in `ScatterPlot.isRowDrawable` /
    `_ScatterPlotSmartLabelsFeature.includeCurrentRowCheck`). This
    scenario covers analyze-then-result; the save-with-datasync →
    reopen → ScatterPlot round-trip invariant is delegated to
    `bio-grok-19928-spec.ts` (chain
    `bug_focused_candidates[GROK-19928]` spans `analyze.md:Step 1` +
    `sequence-activity-cliffs.md:Step 2` +
    `sequence-space.md:Step 2`). Atlas cross-cutting:
    `bio.x.bio-analysis-in-datasync-projects`.
- **Cross-scenario context:** `analyze.md` exercises Sequence Space
  alongside Activity Cliffs and Composition as an umbrella runner;
  `sequence-activity-cliffs.md` is the parallel deeper scenario for
  the Activity-Cliffs top menu. Together `sequence-space.md` +
  `sequence-activity-cliffs.md` realize the two atlas critical
  paths under the analyze-search axis (`bio.cp.sequence-space` and
  `bio.cp.activity-cliffs`).
- **Unresolved ambiguities** (carried in frontmatter
  `unresolved_ambiguities`):
  - The source text says "Change the Similarity metric and Method
    name parameters arbitrarily" — atlas does not enumerate a
    canonical edit set for the `bio.analyze.sequence-space.editor`
    Similarity metric / Method name inputs. The prior run log
    (`sequence-space-run.md`) used UMAP → t-SNE and Hamming →
    Levenshtein (Needleman-Wunsch on FASTA only) as concrete picks,
    but neither atlas nor an operator decision binds those as the
    canonical edit set.
  - The dialog's Macromolecule-column auto-selection contract is
    not declared in atlas (the prior run log notes the dialog
    auto-picks the first Macromolecule column — `Sequence` / `HELM`
    / `MSA` on the three fixtures; on a multi-Macromolecule-column
    dataframe the selection contract is unspecified). Operator
    clarification is needed before Automator can encode a typed
    assertion on the Sequence-Space-editor column selection.
- **Deferred via SCOPE_REDUCTION** (see frontmatter
  `scope_reductions`):
  - SR-01 — the correctness assertion on the edited-parameter run
    (source Step 6) is deferred until atlas or operator supplies a
    concrete Similarity metric / Method name edit set. The
    edit-then-run flow itself is preserved as a step (dialog
    re-opens, inputs are editable, a second embedding ScatterPlot
    docks).
- **`derived_from:` provenance:** atlas entries cited above derive
  from code anchors only — `bio.analyze.sequence-space.*` from
  `package.ts#L237-L789`, `bio.cp.sequence-space` from
  `package.ts#L740`, `bio.cp.fasta-import-via-multiple-entry-paths`
  from `bug-library/bio.yaml#GROK-18616`. No help-doc derivations
  were used (per the binding sourcing rule).

---
{
  "order": 9
}
