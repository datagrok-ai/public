---
feature: bio
sub_features_covered:
  - bio.calculate.identity
  - bio.calculate.similarity
  - bio.calculate.seq-identity
  - bio.calculate.get-region.api
  - bio.analyze.alignment-pairwise
  - bio.calculate.get-region
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-calculate-scoring-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T08:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T08:13:00Z
    spec_runs:
      - spec: bio-calculate-scoring-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 124
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T09:00:00Z
    failure_keys: []
    review_round: 1
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
---

# Bio | Calculate — Identity / Similarity scoring + region API + pairwise alignment

Gate F SR extension scenario (cycle 2026-06-01-bio-migrate-02) that
realizes atlas critical path `bio.cp.identity-scoring` (`priority: p1`)
and extends coverage of the Bio Calculate family beyond the
top-menu-only entries already realized in earlier section scenarios
(`bio.calculate.get-region` and `bio.calculate.get-region.top-menu`
are covered upstream via `bio.flow.get-region`-helper usage). Targets
five net-new bio sub-features on the live covered union and one
already-covered umbrella id (`bio.calculate.get-region`) carried as
the API anchor and for the F-STRUCT-INTERACTION-01 density floor.

The scenario exercises four entry points into the Bio score / region
/ pairwise-alignment surface:

(a) `Bio | Calculate | Identity...` top-menu (atlas
    `bio.calculate.identity`, `package.ts#L1321`) — adds a column of
    fraction-of-matching-monomers vs a reference sequence.
(b) `Bio | Calculate | Similarity...` top-menu (atlas
    `bio.calculate.similarity`, `package.ts#L1336`) — adds a column
    of summed monomer-fingerprint similarities vs a reference.
(c) The public API wrappers `seqIdentity(seq, ref)` (atlas
    `bio.calculate.seq-identity`, `package.ts#L1640`),
    `getRegion(sequence, start, end, name)` (atlas
    `bio.calculate.get-region.api`, `package.ts#L473`), and
    `sequenceAlignment(alignType, alignTable, gap, seq1, seq2)`
    (atlas `bio.analyze.alignment-pairwise`, `package.ts#L435`)
    invoked via `grok.functions.call('Bio:<name>', ...)`.
(d) The umbrella parent `bio.calculate.get-region` anchor —
    exercised at API level (the .api child is one entry point;
    the .top-menu sibling is covered by `composition-analysis.md` /
    `convert.md` upstream).

Per atlas (`bio.cp.identity-scoring`, `derived_from:
package.ts#L1321`, `priority: p1`), the canonical
identity/similarity-scoring path is a `coverage_type: regression`
surface — the score-column outputs feed Sequence Space /
Activity Cliffs analyze surfaces and any change to the
`calculateScoresWithEmptyValues` empty-value handling can regress
multiple downstream consumers without a direct UI symptom.

## Setup

Datasets (canonical Bio test fixtures):

- `System.AppData/Bio/tests/filter_HELM.csv` — HELM notation
  (detector sets `units=helm`); used by Scenarios 1–2 and 4. HELM
  peptides expose the broadest monomer-library lookup surface, so
  Identity / Similarity scoring exercises the
  `calculateScoresWithEmptyValues` wrapper against real-world
  monomer-fingerprint paths.
- `System.AppData/Bio/tests/filter_FASTA.csv` — FASTA-notation
  peptide sequences; used by Scenario 3 to exercise `getRegion` API
  on a FASTA-notation Macromolecule column (`bio.calculate.get-region.api`
  works across notations via `ISeqHandler.getRegion()`).
- `System.AppData/Bio/tests/filter_FASTA.csv` is reused by
  Scenario 5 for the `sequenceAlignment` pairwise alignment — two
  short FASTA sequences picked as `seq1`/`seq2` parameters.

`grok.functions.call('Bio:seqIdentity', {seq, ref})` and
`grok.functions.call('Bio:getRegion', {sequence, start, end, name})`
and `grok.functions.call('Bio:sequenceAlignment', {alignType,
alignTable, gap, seq1, seq2})` are the public API entry points;
their `@grok.decorators.func` registrations live alongside the
top-menu functions in `public/packages/Bio/src/package.ts`.

Helper handle: the `bio.flow.get-region` helper-registry entry
provides the dataset-open + region-extraction skeleton shared with
upstream scenarios; this scenario reuses the dataset-open prelude
but invokes the API directly rather than the top-menu Extract
Region dialog.

## Scenarios

### Scenario 1 — Identity scoring via top-menu

Steps:

1. Open `System.AppData/Bio/tests/filter_HELM.csv` from the
   Files browser. The Macromolecule detector classifies the
   sequence column and sets `units=helm`.
2. From the ribbon, run
   **Bio > Calculate > Identity...** (atlas
   `bio.calculate.identity`, `package.ts#L1321`). The dialog
   opens prefilled with the active table and Macromolecule
   column.
3. In the dialog `Reference` field, paste the first row's
   sequence (the canonical Scoring-tests reference shape
   `PEPTIDE1{...}$$$$`, see `Bio/src/tests/scoring.ts#L23`).
   Click **OK**.
4. The transform produces a new score column appended to the
   table.

Expected:

- After step 3 the dialog closes cleanly and the new score
  column appears on the dataframe.
- The first row's score is `1.0` (matching reference == self
  identity by definition); subsequent rows carry scores in the
  closed interval `[0.0, 1.0]` corresponding to the
  fraction-of-matching-monomers metric described at
  `package.ts#L1322` `description: 'Adds a column with fraction
  of matching monomers'`.
- No error balloon appears (`isErrorBallon` returns false).

### Scenario 2 — Similarity scoring via top-menu (HELM monomer fingerprints)

Steps:

1. On the same `filter_HELM.csv` table view from Scenario 1
   (or re-open it if Scenario 1 was run independently), run
   **Bio > Calculate > Similarity...** (atlas
   `bio.calculate.similarity`, `package.ts#L1336`). The dialog
   opens prefilled.
2. In the dialog `Reference` field, paste the first row's
   sequence. Click **OK**.
3. The transform produces a new similarity score column
   appended to the table.

Expected:

- After step 2 the dialog closes cleanly and the new
  similarity column appears alongside any identity column
  from Scenario 1.
- The first row's similarity score is `1.0` (self-similarity);
  subsequent rows carry monotonically meaningful float scores
  produced by the summed monomer-fingerprint similarity metric
  documented at `package.ts#L1338` `description: 'Adds a
  column with similarity scores, calculated as sum of monomer
  fingerprint similarities'`.
- No error balloon appears.

### Scenario 3 — Region extraction via `getRegion` API call (Macromolecule column, FASTA notation)

Steps:

1. Open `System.AppData/Bio/tests/filter_FASTA.csv` from the
   Files browser. The Macromolecule detector sets `units=fasta`
   on the sequence column.
2. Invoke the `getRegion` API directly via
   `await grok.functions.call('Bio:getRegion', {sequence:
   <col>, start: 0, end: 4, name: 'region_0_4'})`. The function
   is the public API surfaced from `@grok.decorators.func({name:
   'getRegion', ...})` at `package.ts#L473`; it delegates to
   `ISeqHandler.getRegion()` (atlas
   `bio.calculate.get-region`, parent anchor) and returns the
   positional sub-region as a new Macromolecule column.
3. Append the returned column to the active dataframe.

Expected:

- The `getRegion` call returns a `DG.Column<string>` whose cells
  contain the first four monomers of each source sequence,
  preserving the FASTA notation tags (`units=fasta`, semType =
  `Macromolecule`).
- The returned column name is `region_0_4` (the `name`
  parameter).
- The umbrella `bio.calculate.get-region` parent surface is
  exercised through `ISeqHandler.getRegion()` (atlas anchor
  `package.ts#L473`); the `.api` child registration is the entry
  point under test here.
- No error balloon appears.

### Scenario 4 — `seqIdentity` single-pair API call (empty-input contract)

Steps:

1. With `filter_HELM.csv` open (or any dataframe — `seqIdentity`
   is a single-pair function, not a column transform), call
   `await grok.functions.call('Bio:seqIdentity', {seq: <first
   row's HELM sequence>, ref: <same first row's HELM
   sequence>})`. Atlas `bio.calculate.seq-identity`,
   `package.ts#L1640`.
2. Call `await grok.functions.call('Bio:seqIdentity', {seq:
   '', ref: <first row's HELM sequence>})` to exercise the
   empty-input branch documented on the atlas
   (`returns null when seq is empty`).
3. Call `await grok.functions.call('Bio:seqIdentity', {seq:
   <row 1 sequence>, ref: <row 0 sequence>})` for a cross-row
   pair to confirm the function returns a non-trivial float.

Expected:

- Step 1 returns `1.0` (self-identity).
- Step 2 returns `null` per the atlas contract — empty `seq`
  short-circuits to null rather than throwing, and downstream
  code paths must tolerate the null cell value (this is the
  empty-value branch handled in
  `utils/calculate-scores.ts :: calculateScoresWithEmptyValues`).
- Step 3 returns a numeric float in `[0.0, 1.0]`.
- No error balloon appears across the three invocations.

### Scenario 5 — Pairwise alignment via `sequenceAlignment` API (Needleman-Wunsch)

Steps:

1. With `filter_FASTA.csv` open (re-opened from Scenario 3 if
   needed), pick two short FASTA-notation sequences from the
   Macromolecule column (rows 0 and 1).
2. Invoke
   `await grok.functions.call('Bio:sequenceAlignment',
   {alignType: 'global', alignTable: 'BLOSUM62', gap: -10,
   seq1: <row 0 sequence>, seq2: <row 1 sequence>})` (atlas
   `bio.analyze.alignment-pairwise`, `package.ts#L435`,
   wrapping the `SequenceAlignment` Needleman-Wunsch /
   Smith-Waterman engine from `src/seq_align.ts`).
3. Repeat with `alignType: 'local'` (Smith-Waterman) and any
   other BLOSUM matrix — e.g. `BLOSUM45` — to exercise the
   matrix-selection branch on the `SequenceAlignment`
   constructor.

Expected:

- Both invocations return a non-null aligned-sequences result
  (the `SequenceAlignment` engine emits a pair of gap-padded
  aligned sequence strings; the precise shape is implementation-
  defined and the test asserts non-null, non-empty, and length
  >= max(len(seq1), len(seq2))).
- The `global` (Needleman-Wunsch) and `local` (Smith-Waterman)
  branches both run without throwing; the BLOSUM45 / BLOSUM62
  matrix-name branch on the alignment-engine constructor is
  exercised across the two invocations.
- No error balloon appears.

## Notes

- atlas entry derived from sub_features `bio.calculate.{identity,
  similarity, seq-identity, get-region.api}` and
  `bio.analyze.alignment-pairwise`; the top-menu steps in
  Scenarios 1–2 derive from atlas critical_path
  `bio.cp.identity-scoring` (`priority: p1`,
  `derived_from: public/packages/Bio/src/package.ts#L1321`).
- target_layer rationale: `playwright` — Scenarios 1 and 2 are
  ribbon-driven top-menu UI dialogs (`Bio > Calculate > Identity`
  and `Bio > Calculate > Similarity`); the score-column output
  must be observable on the painted grid for the dialog-close
  contract. Scenarios 3–5 are API-call shape inside the same
  Playwright session via `evaluate` against
  `grok.functions.call(...)` — combining them into one scenario
  file preserves the section-wide F-STRUCT-DENSITY-01 average
  (6 sub-features over the scenario file) and the
  F-STRUCT-INTERACTION-01 floor (six sub-features in
  interaction). A pure-apitest factoring is possible for the
  API-only branches but would split the Calculate family across
  two scenario files for no additional value.
- coverage_type rationale: `regression` — score-column outputs
  feed downstream analyze surfaces (Sequence Space distance
  inputs, Activity Cliffs activity-pair detection, MSA
  per-cluster scoring) and the
  `calculateScoresWithEmptyValues` empty-value path was
  observed to regress on prior cycles (no specific bug filed
  but covered by the regression-class contract). Atlas
  critical_path `bio.cp.identity-scoring` carries `priority:
  p1` which maps to `coverage_type: regression` per the
  skill's STEP E heuristic.
- Sub-features covered:
  - `bio.calculate.identity` (`package.ts#L1321`) — Identity
    top-menu transform; covered by Scenario 1.
  - `bio.calculate.similarity` (`package.ts#L1336`) —
    Similarity top-menu transform; covered by Scenario 2.
  - `bio.calculate.seq-identity` (`package.ts#L1640`) —
    single-pair `seqIdentity` API; covered by Scenario 4
    (including the empty-input null-return contract).
  - `bio.calculate.get-region.api` (`package.ts#L473`) —
    public `getRegion` API; covered by Scenario 3.
  - `bio.analyze.alignment-pairwise` (`package.ts#L435`) —
    `sequenceAlignment` pairwise engine; covered by
    Scenario 5.
  - `bio.calculate.get-region` (`package.ts#L473`) — parent
    anchor for the region surface; exercised through the
    `getRegion` API call in Scenario 3 (the `.api` child is
    one entry point into the umbrella).
- Net-new ids vs `live_covered_union` (this cycle's
  authoritative covered set of 56 ids):
  - `bio.calculate.identity` — net-new.
  - `bio.calculate.similarity` — net-new.
  - `bio.calculate.seq-identity` — net-new.
  - `bio.calculate.get-region.api` — net-new.
  - `bio.analyze.alignment-pairwise` — net-new.
  - `bio.calculate.get-region` — already in
    `live_covered_union` (anchor, not net-new).
  - net_new = 5 ids; satisfies the SR-loop progress-sensitive
    bound (delta > 0) and the STEP C net-new refusal.
    Projected coverage after this iteration's merge: 59/99
    (~59.6%), advancing toward the 70% threshold while the
    remaining first-pass scenarios are authored on subsequent
    iterations.
- Manual-only subset: none of the six covered sub_features
  appear in atlas `manual_only[]` (verified against atlas rev
  3 `manual_only[]` list — none of `bio.calculate.*` or
  `bio.analyze.alignment-pairwise` is flagged manual_only).
- Deferrals: none. All five scenarios are observable inside a
  single Playwright session (UI ribbon for Scenarios 1–2;
  `grok.functions.call` via `evaluate` for Scenarios 3–5).
- Bug-context (`related_bugs`): none. No atlas
  `known_issues[]` / `edge_cases[]` entry maps to the
  Calculate scoring family directly; the GROK-12164 dispatch
  regression is renderer-scope, not score-scope. `related_bugs:
  []` per atlas.
- The Bio section has no grok-browser ref-doc verb-form H2
  matching the citation regex (see
  `references/bio.md` headings audit by F-UI-COVERAGE-01 for
  this cycle), so `## Notes` citation pointers reference atlas
  entries only — no `See: <path>#<heading>` citation form
  applies here.
- This scenario covers 6 sub_features (`F-STRUCT-DENSITY-01`
  floor: 2; `F-STRUCT-INTERACTION-01` floor: 3 in a
  multi-sub_feature scenario — satisfied).

---
{
  "order": 22
}
