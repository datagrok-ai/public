---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [bio.cp.identity-scoring]
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
---

# Bio | Calculate — Identity / Similarity scoring + region API + pairwise alignment

Covers the Bio scoring and alignment surface: the **Identity...** and
**Similarity...** top-menu actions (each adds a score column comparing
every sequence in a column against a reference sequence), plus the
underlying public API functions — `seqIdentity`, `getRegion`, and
`sequenceAlignment` (pairwise Needleman-Wunsch / Smith-Waterman
alignment). These scoring outputs feed the Sequence Space and Activity
Cliffs analyses, so a regression here can silently break those
downstream features without an obvious UI symptom.

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
   **Bio > Calculate > Identity...**. The dialog
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
   **Bio > Calculate > Similarity...**. The dialog
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
   `ISeqHandler.getRegion()` and returns the
   positional sub-region as a new Macromolecule column.
3. Append the returned column to the active dataframe.

Expected:

- The `getRegion` call returns a `DG.Column<string>` whose cells
  contain the first four monomers of each source sequence,
  preserving the FASTA notation tags (`units=fasta`, semType =
  `Macromolecule`).
- The returned column name is `region_0_4` (the `name`
  parameter).
- The underlying `ISeqHandler.getRegion()` region-extraction
  surface is exercised through this API call.
- No error balloon appears.

### Scenario 4 — `seqIdentity` single-pair API call (empty-input contract)

Steps:

1. With `filter_HELM.csv` open (or any dataframe — `seqIdentity`
   is a single-pair function, not a column transform), call
   `await grok.functions.call('Bio:seqIdentity', {seq: <first
   row's HELM sequence>, ref: <same first row's HELM
   sequence>})`.
2. Call `await grok.functions.call('Bio:seqIdentity', {seq:
   '', ref: <first row's HELM sequence>})` to exercise the
   empty-input branch (`returns null when seq is empty`).
3. Call `await grok.functions.call('Bio:seqIdentity', {seq:
   <row 1 sequence>, ref: <row 0 sequence>})` for a cross-row
   pair to confirm the function returns a non-trivial float.

Expected:

- Step 1 returns `1.0` (self-identity).
- Step 2 returns `null` — empty `seq`
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
   seq1: <row 0 sequence>, seq2: <row 1 sequence>})` (wrapping
   the `SequenceAlignment` Needleman-Wunsch / Smith-Waterman
   engine from `src/seq_align.ts`).
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

---
{
  "order": 22
}
