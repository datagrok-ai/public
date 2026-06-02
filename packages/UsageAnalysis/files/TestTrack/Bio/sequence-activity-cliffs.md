---
feature: bio
sub_features_covered:
  - bio.analyze.activity-cliffs
  - bio.analyze.activity-cliffs.top-menu
  - bio.analyze.activity-cliffs.editor
  - bio.analyze.activity-cliffs.transform
  - bio.analyze.activity-cliffs.init
target_layer: playwright
coverage_type: regression
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
    scope_reduction_proposal: |
      Adjudicate SR-01 (check: A-CONT-01) as SCOPE_REDUCTION. The
      scenario body has already applied the reduction in place: source
      Step 5's non-deterministic "Change the Similarity and Method name
      parameters arbitrarily" is preserved as the edit-then-run flow in
      Step 6 (dialog must re-open, inputs must accept new selections via
      the select widget, exercising the editor input re-binding
      contract), while Step 7's correctness assertion on the
      edited-parameter SALI distribution is explicitly deferred. The
      observable assertion that survives is structural — a second
      "Activity cliffs" ScatterPlot must dock, distinct from the first
      run (total viewer count = 2). The cited dependency is real and
      atlas-level: atlas sub-feature bio.analyze.activity-cliffs.editor
      (L380, SeqActivityCliffsEditor at package.ts#L268) describes the
      dialog inputs (Similarity metric, Method, similarity cutoff) but
      does NOT enumerate which Similarity / Method combinations form a
      canonical edit set, and the unresolved_ambiguities[] field
      acknowledges that the prior run log used Hamming->Levenshtein and
      UMAP->t-SNE as concrete picks but neither atlas nor an operator
      decision binds those as canonical. Bumping SR-01.verdict_status
      from null to SCOPE_REDUCTION; once atlas or operator supplies a
      concrete edit set, a follow-up scenario revision can promote the
      deferred assertion into a typed SALI-distribution check.
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

Integration scenario for the `Bio | Analyze | Activity Cliffs...` top
menu (atlas `bio.analyze.activity-cliffs.top-menu` — `package.ts#L537`)
across three canonical Macromolecule notations (FASTA, HELM, MSA), with
two run modes (default parameters, then edited Similarity + Method
name). Verifies the multi-subsystem path: editor dialog
(`bio.analyze.activity-cliffs.editor`) → embedding compute
(`bio.analyze.activity-cliffs.transform`) → ScatterPlot viewer with
cliff overlay (`bio.analyze.activity-cliffs.init`).

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
assertion is deferred (see SR-01).

## Scenarios

### Scenario 1 — Open dataset and dispatch Activity Cliffs with defaults

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector (atlas `bio.detector`) classifies the sequence column
   synchronously; the table view opens.

2. On the menu ribbon, open **Bio** > **Analyze** > **Activity
   Cliffs...**. The Activity Cliffs dialog opens (atlas
   `bio.analyze.activity-cliffs.editor`,
   `SeqActivityCliffsEditor` — `package.ts#L268`). The dialog
   collects: Macromolecule column, Activities (numeric), Similarity
   metric, Method, similarity cutoff.

3. Click **OK** to run with default parameters. The transform
   dispatches (atlas `bio.analyze.activity-cliffs.transform`,
   `seqActivityCliffsTransform` — `package.ts#L658`), computing
   embeddings and stamping the `seqActivityCliffsParams` tag on the
   table.

4. Verify the result viewer is added:
   - A ScatterPlot viewer titled "Activity cliffs" docks (atlas
     `bio.analyze.activity-cliffs.init`,
     `seqActivityCliffsInitFunction` — `package.ts#L625`) with the
     cliff overlay drawn on top of the sequence-space embedding.
   - Embedding columns (X / Y) and the SALI scoring column are
     appended to the source dataframe.

### Scenario 2 — Re-open dialog and run with edited Similarity / Method

For each dataset (same matrix as Scenario 1):

5. On the menu ribbon, re-open **Bio** > **Analyze** > **Activity
   Cliffs...** The dialog opens with the prior run's defaults.

6. In the dialog, change the **Similarity** metric and **Method name**
   inputs (the source text says "arbitrarily" — atlas does not pin a
   canonical edit set, see SR-01 and `unresolved_ambiguities`). The
   inputs must remain editable (the dialog accepts new selections via
   the `<select>` widget surface) — this exercises the editor input
   re-binding contract.

7. Click **OK** to run with the edited parameters. The transform
   re-dispatches; verify a second "Activity cliffs" ScatterPlot
   docks, distinct from the first run's viewer (total Activity Cliffs
   viewer count = 2). The edited-parameter result correctness
   assertion (e.g. that a `Levenshtein + t-SNE` run produces a
   meaningfully different SALI distribution than a `Hamming + UMAP`
   run) is deferred — see SR-01.

## Notes

- **Source-text fixes silently applied** during migration (recorded in
  frontmatter `source_text_fixes`):
  - The original scenario placed Activity Cliffs under **Bio > Search
    > Sequence Activity Cliffs**. Atlas registers this surface as
    `bio.analyze.activity-cliffs.top-menu` under **Bio | Analyze |
    Activity Cliffs...** (`package.ts#L537`); the prior run log
    (`sequence-activity-cliffs-run.md`) confirms the menu path "Bio >
    Search" does not exist on the live platform and that the live
    function lives at "Bio > Analyze > Activity Cliffs...". The chain
    flagged this as an unresolved ambiguity for operator confirmation;
    atlas + run-log evidence converge, so the migrated scenario uses
    the atlas-canonical path.
  - The original scenario referred to the function as "Sequence
    Activity Cliffs" (verbatim function name, no ellipsis). The
    atlas-canonical menu label is "Activity Cliffs..." with the
    trailing ellipsis denoting a dialog-opening function; the function
    is internally registered as "Sequence Activity Cliffs" but its
    user-facing menu label drops "Sequence" and adds the ellipsis. The
    migrated scenario uses the user-facing label "Activity Cliffs..."
    consistent with atlas `bio.analyze.activity-cliffs.top-menu`
    (`package.ts#L537`).
  - The implicit 3-dataset matrix was listed in Setup but not
    explicitly framed; surfaced here under Setup as a 3 × 2 matrix.
  - The original scenario named the datasets `sample_FASTA.csv`,
    `sample_HELM.csv`, `sample_MSA.csv` (with the `sample_` prefix and
    no directory path). The repo's Bio package ships no such files —
    the closest `sample_`-prefixed file in `Bio/files/` is
    `sample_MSA_data.csv` (tests/), which is unrelated. The canonical
    Activity-Cliffs fixtures are `filter_FASTA.csv`, `filter_HELM.csv`,
    `filter_MSA.csv` under `Bio/files/tests/` (the `filter_` family is
    purpose-built for Activity-Cliffs / Sequence-Space and is what the
    sister `analyze-spec.ts`, `sequence-space-spec.ts`, and
    `sequence-activity-cliffs-spec.ts` resolve to via
    `System.AppData/Bio/tests/`). Migrated scenario pins the dataset
    family to `filter_*.csv` to match the existing spec corpus and to
    keep the 3-dataset matrix exercising notation-agnostic behaviour
    (FASTA / HELM / MSA). The prior run log
    (`sequence-activity-cliffs-run.md`) used `samples/FASTA.csv` etc.
    — a different family from the migrated choice; the divergence
    between run-log empirical resolution and migrated fixture choice
    is acknowledged here. Operator clarification welcomed if the
    `samples/` family is preferred over the purpose-built `tests/`
    family.
  - The original scenario gave only bare filenames with no directory
    path. The migrated scenario pins the path to
    `System.AppData/Bio/tests/` to match the `filter_*.csv` family
    selection above and to align with the sister specs' resolution
    pattern.
- **Sub-features covered:**
  - `bio.analyze.activity-cliffs` (atlas L346) — feature root.
  - `bio.analyze.activity-cliffs.top-menu` (atlas L355) — `Bio |
    Analyze | Activity Cliffs...` dispatch (`package.ts#L537`).
  - `bio.analyze.activity-cliffs.editor` (atlas L380) —
    `SeqActivityCliffsEditor` dialog (`package.ts#L268`).
  - `bio.analyze.activity-cliffs.transform` (atlas L364) —
    `seqActivityCliffsTransform` embedding compute
    (`package.ts#L658`).
  - `bio.analyze.activity-cliffs.init` (atlas L372) —
    `seqActivityCliffsInitFunction` cliff-overlay paint on the
    ScatterPlot viewer (`package.ts#L625`).
  Maps directly onto atlas critical path `bio.cp.activity-cliffs`
  (p0, `derived_from: public/packages/Bio/src/package.ts#L537`).
- **Related bugs** surfaced for downstream awareness; per-bug
  cross-cutting invariants are delegated to chain-level
  `bug_focused_candidates`:
  - GROK-19150 — Activity Cliffs multi-instance state-isolation gap
    (clicks on the second analysis route to the first analysis's
    `cliffs_resulting_grid`). This scenario covers single-instance
    Activity Cliffs only; the click-routing isolation invariant
    requires opening ≥ 2 analyses with distinct methods and is
    delegated to `bio-grok-19150-spec.ts` (chain
    `bug_focused_candidates[GROK-19150]` spans
    `analyze.md:Step 1` + `sequence-activity-cliffs.md:Step 2`).
    Atlas cross-cutting:
    `bio.x.activity-cliffs-multi-instance-isolation`.
  - GROK-19928 — Sequence Space and Activity Cliffs analysis output
    not properly saved/restored in datasync projects (reopen crashes
    with `NullError` in `ScatterPlot.isRowDrawable` /
    `_ScatterPlotSmartLabelsFeature.includeCurrentRowCheck`). This
    scenario covers analyze-then-result; the save-with-datasync →
    reopen → ScatterPlot round-trip invariant is delegated to
    `bio-grok-19928-spec.ts` (chain `bug_focused_candidates[GROK-19928]`
    spans `analyze.md:Step 1` + `sequence-activity-cliffs.md:Step 2`
    + `sequence-space.md:Step 2`). Atlas cross-cutting:
    `bio.x.bio-analysis-in-datasync-projects`.
  - GROK-16111 — Bio analyze viewers operating on the current row
    (Sequence Similarity Search, Diversity Search, Activity Cliffs)
    must reject empty / null input via balloon. This scenario uses
    non-empty test fixtures; the empty-input rejection invariant is
    an atlas-cross-cutting edge case
    (`bio.x.empty-input-on-row-viewers`) covered separately and is
    not asserted here.
- **Cross-scenario context:** `analyze.md` exercises Activity Cliffs
  alongside Sequence Space and Composition as an umbrella runner;
  `sequence-space.md` is the parallel deeper scenario for the
  Sequence-Space top menu. Together
  `sequence-activity-cliffs.md` + `sequence-space.md` realize the
  two atlas critical paths under the analyze-search axis
  (`bio.cp.activity-cliffs` and `bio.cp.sequence-space`).
- **Unresolved ambiguities** (carried in frontmatter
  `unresolved_ambiguities`):
  - The source text says "Change the Similarity and Method name
    parameters arbitrarily" — atlas does not enumerate a canonical
    edit set for the `bio.analyze.activity-cliffs.editor` Similarity
    / Method inputs. The prior run log
    (`sequence-activity-cliffs-run.md`) used Hamming → Levenshtein
    and UMAP → t-SNE as concrete picks, but neither atlas nor an
    operator decision binds those as the canonical edit set.
  - The dialog's **Activities** input defaults to the first numeric
    column (per the prior run log, the `Length` column was picked
    rather than the `Activity` column on the test fixtures). The
    default-column-detection contract for Activity Cliffs is not
    declared in atlas; operator clarification is needed before
    Automator can encode a typed assertion on the Activities-column
    selection.
- **Deferred via SCOPE_REDUCTION** (see frontmatter
  `scope_reductions`):
  - SR-01 — the correctness assertion on the edited-parameter run
    (source Step 6) is deferred until atlas or operator supplies a
    concrete Similarity / Method edit set. The edit-then-run flow
    itself is preserved as a step (dialog re-opens, inputs are
    editable, a second cliff result docks).
- **`derived_from:` provenance:** atlas entries cited above derive
  from code anchors only — `bio.analyze.activity-cliffs.*` from
  `package.ts#L268-L658`, `bio.cp.activity-cliffs` from
  `package.ts#L537`. No help-doc derivations were used (per the
  binding sourcing rule).

---
{
  "order": 8
}
