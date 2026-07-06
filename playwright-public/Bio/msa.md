---
feature: bio
sub_features_covered:
  - bio.analyze.msa
  - bio.analyze.msa.dialog
  - bio.analyze.msa.align-sequences
  - bio.rendering.column-header
  - bio.detector
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/msa.md
migration_date: 2026-05-31
source_text_fixes:
  - step-5-trailing-whitespace-trimmed
  - step-6-assertion-tail-expanded-into-three-verification-steps
candidate_helpers: []
realized_as:
  - msa-spec.ts
unresolved_ambiguities:
  - cluster-column-input-type-contract-integer-vs-categorical
scope_reductions: []
related_bugs:
  - GROK-18474
  - GROK-15176
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T09:30:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T10:00:00Z
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
    timestamp: 2026-06-01T14:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T23:50:00Z
    spec_runs:
      - spec: msa-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 53
        failure_keys: []
---

# Bio | Analyze | MSA — canonical FASTA (kalign) integration

Integration scenario for the canonical Multiple Sequence Alignment path
(atlas `bio.analyze.msa`, `bio.analyze.msa.dialog`,
`bio.analyze.msa.align-sequences` — `package.ts#L968`). Runs kalign WASM
on a canonical FASTA peptide dataset with a per-row Cluster column, then
verifies the result aligned column is added, the MSA column-header
renderer is wired (atlas `bio.rendering.column-header`), and the
per-cluster alignment invariant holds.

The non-canonical (HELM → PepSeA Docker engine) counterpart lives in
`pepsea.md`; this scenario covers the kalign WASM branch on FASTA.

## Setup

Dataset: `System.AppData/Bio/tests/filter_FASTA.csv` (canonical FASTA
peptide fixture per chain analyzer for this section). The dataset is
opened so the Macromolecule detector (atlas `bio.detector`) classifies
the sequence column synchronously — the `Bio | Analyze | MSA...` top
menu is only reachable when a Macromolecule column is present on the
active table view.

A second column is then added via formula `RandBetween(0, 5)` to serve
as the per-row Cluster input to the MSA dialog. The chain analyzer
flagged the integer-vs-categorical contract of the dialog's Cluster
input as an unresolved ambiguity (carried in frontmatter
`unresolved_ambiguities`); this scenario follows the source text and
supplies the formula-produced integer column directly.

## Scenarios

### Scenario 1 — Open dataset, add cluster column, dispatch MSA dialog

1. Open `System.AppData/Bio/tests/filter_FASTA.csv`. The Macromolecule
   detector classifies the sequence column (atlas `bio.detector`); the
   table view opens.

2. Add a new column with formula `RandBetween(0, 5)`. This column will
   serve as the per-row Cluster input in step 4.

3. On the menu ribbon, open **Bio** > **Analyze** > **MSA...**. The MSA
   dialog opens (atlas `bio.analyze.msa.dialog`,
   `multipleSequenceAlignmentDialog` — `package.ts#L968`).

4. In the MSA dialog, set the **Cluster** input to the column produced
   in step 2 (the `RandBetween(0, 5)` formula column).

### Scenario 2 — Alignment parameters button toggles dialog inputs

5. Click the **Alignment parameters** button in the MSA dialog. Verify
   the button adds input parameters to the dialog (the dialog grows the
   per-engine alignment parameters surface). Click again to confirm the
   button toggles the parameter surface; the dialog returns to its
   prior shape.

### Scenario 3 — Run MSA and verify result column

6. Click **OK** to run alignment. MSA dispatches by alphabet: canonical
   (Peptide) routes to the kalign WASM engine (atlas
   `bio.analyze.msa.align-sequences` — `package.ts#L990`).

7. Verify the result aligned column is added to the table (a new MSA
   column appears alongside the source sequence column).

8. Verify the MSA column-header renderer is set on the result column
   (atlas `bio.rendering.column-header` — WebLogo header + conservation
   tracks paint on the aligned column; the renderer dispatches on the
   column's MSA flag).

9. Verify the per-cluster alignment invariant holds: for each Cluster
   value in the `RandBetween(0, 5)` column, the aligned sequences in
   that cluster are of the same length (kalign aligns within each
   Cluster group; cross-cluster lengths may differ).

## Notes

- **Source-text fixes silently applied** during migration (recorded in
  frontmatter `source_text_fixes`):
  - Step 5 had a trailing run of spaces after "properly"; trimmed.
  - Original step 6 said only "Check the new column." — vague. The
    assertion tail line ("Everything is good if the new MSA column is
    added, the renderer for the column is set, and the sequences within
    a cluster are of the same length.") has been promoted into three
    explicit verification steps (Scenario 3, steps 7-9 above) so the
    assertion surface is decidable. Per D-STEP-02 the expected-result
    line is preserved verbatim in semantic content; no assertion was
    dropped.
- **Sub-features covered:**
  - `bio.analyze.msa` (atlas L422) — feature root.
  - `bio.analyze.msa.dialog` (atlas L431) — `Bio | Analyze | MSA...`
    top-menu surface (`multipleSequenceAlignmentDialog`).
  - `bio.analyze.msa.align-sequences` (atlas L441) — canonical kalign
    engine path (FASTA peptides are canonical → kalign WASM).
  - `bio.rendering.column-header` (atlas — MSA-aware column-header
    renderer) — exercised in step 8 verification.
  - `bio.detector` (atlas L122) — touched in setup (Macromolecule
    classification of the FASTA sequence column on dataset open).
  Maps directly onto atlas critical path `bio.cp.msa-canonical` (p0,
  `derived_from: public/packages/Bio/src/package.ts#L968`).
- **Related bugs** surfaced for downstream awareness; per-bug
  cross-cutting coverage is delegated to chain-level
  `bug_focused_candidates`:
  - GROK-18474 — MSA column-header click handler crashes on FASTA-
    aligned data (affects `bio.rendering.column-header`, exercised by
    step 8). The dock-on-click invariant is the cross-notation regression
    check (atlas `bio.x.msa-header-click-cross-notation`); this
    scenario asserts renderer presence, not click handler behavior.
  - GROK-15176 — Bio's to-atomic-level produces molfiles with invalid
    isotope on heavy atoms (affects `bio.analyze.msa` +
    `bio.analyze.msa.dialog` as upstream producers). Chain proposes
    `bio-grok-15176-spec.ts` spanning `msa.md:Step 3`, `pepsea.md:Step
    3`, `convert.md:Step 2`; cross-cutting MSA → To Atomic Level
    invariant is not asserted here. See chain
    `bug_focused_candidates[GROK-15176]`.
- **Helper coverage:** `bio.flow.msa` already registered in
  `helpers-registry.yaml:219` →
  `.claude/skills/grok-browser/references/bio.md:139`. The single
  registered helper covers both `msa.md` (this scenario) and
  `pepsea.md` (non-canonical counterpart). No new candidate helpers
  proposed.
- **Counterpart scenarios:** `pepsea.md` covers the non-canonical
  (HELM → PepSeA Docker engine) branch on the same MSA top-menu
  surface; atlas `bio.cp.msa-pepsea` (p1) is its canonical critical
  path. Together `msa.md` + `pepsea.md` realize both atlas-declared
  MSA critical paths.
- **Unresolved ambiguity** (carried in frontmatter
  `unresolved_ambiguities`): chain analyzer (`bio.yaml` ambiguity for
  `msa.md:Step 4`) noted the source does not specify whether the MSA
  dialog's Cluster input accepts an integer column (such as the
  `RandBetween(0, 5)` column produced in step 2) directly, or whether
  it expects a categorical/string column. This scenario follows the
  source text and supplies the integer column as-is; operator
  clarification of the Cluster-input column-type contract is needed
  before Automator can encode a typed assertion.
- **`derived_from:` provenance:** the atlas entries cited above were
  derived from code anchors only — `bio.analyze.msa.dialog` and
  `bio.cp.msa-canonical` both derive from `package.ts#L968`,
  `bio.analyze.msa.align-sequences` from `package.ts#L990`. No
  help-doc derivations were used (per the binding sourcing rule).

---
{
  "order": 6
}
