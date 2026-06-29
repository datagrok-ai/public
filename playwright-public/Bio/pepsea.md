---
feature: bio
sub_features_covered:
  - bio.analyze.msa
  - bio.analyze.msa.dialog
  - bio.engines.msa-pepsea
  - bio.rendering.column-header
  - bio.detector
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/pepsea.md
migration_date: 2026-05-31
source_text_fixes:
  - step-5-trailing-whitespace-trimmed
  - step-6-assertion-tail-expanded-into-three-verification-steps
candidate_helpers: []
unresolved_ambiguities:
  - menu-path-bio-msa-vs-bio-analyze-msa
  - cluster-column-input-type-contract-integer-vs-categorical
scope_reductions: []
related_bugs:
  - GROK-15176
realized_as:
  - pepsea-spec.ts
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
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T16:05:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T15:49:00Z
    spec_runs:
      - spec: pepsea-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 85
        failure_keys: []
---

# Bio | Analyze | MSA — non-canonical HELM (PepSeA Docker) integration

Integration scenario for the non-canonical Multiple Sequence Alignment
path (atlas `bio.analyze.msa`, `bio.analyze.msa.dialog`,
`bio.engines.msa-pepsea` — `package.ts#L968`, `#L1003`). Runs the PepSeA
Docker engine on a HELM (non-canonical peptide alphabet) dataset with a
per-row Cluster column, then verifies the result aligned column is
added, the MSA column-header renderer is wired (atlas
`bio.rendering.column-header`), and the per-cluster alignment invariant
holds.

The canonical (FASTA → kalign WASM) counterpart lives in `msa.md`; this
scenario covers the PepSeA Docker-engine branch on HELM. Together
`msa.md` + `pepsea.md` realize both atlas-declared MSA critical paths
(`bio.cp.msa-canonical` p0 and `bio.cp.msa-pepsea` p1).

## Setup

Dataset: `System.AppData/Bio/tests/filter_HELM.csv` (HELM macromolecule
fixture per chain analyzer for this section). The dataset is opened so
the Macromolecule detector (atlas `bio.detector`) classifies the
sequence column synchronously — the `Bio | Analyze | MSA...` top menu
is only reachable when a Macromolecule column is present on the active
table view, and HELM-aligned data routes MSA through the non-canonical
PepSeA Docker engine (atlas `bio.engines.msa-pepsea`,
`package.ts#L1003`) instead of the canonical kalign WASM path.

A second column is then added via formula `RandBetween(0, 5)` to serve
as the per-row Cluster input to the MSA dialog. The chain analyzer
flagged the integer-vs-categorical contract of the dialog's Cluster
input as an unresolved ambiguity (carried in frontmatter
`unresolved_ambiguities`); this scenario follows the source text and
supplies the formula-produced integer column directly.

## Scenarios

### Scenario 1 — Open dataset, add cluster column, dispatch MSA dialog

1. Open `System.AppData/Bio/tests/filter_HELM.csv`. The Macromolecule
   detector classifies the sequence column as HELM (atlas
   `bio.detector`); the table view opens.

2. Add a new column with formula `RandBetween(0, 5)`. This column will
   serve as the per-row Cluster input in step 4.

3. On the menu ribbon, open **Bio** > **Analyze** > **MSA...**. The MSA
   dialog opens (atlas `bio.analyze.msa.dialog`,
   `multipleSequenceAlignmentDialog` — `package.ts#L968`). Source text
   listed the path as **Bio** > **MSA** (no Analyze submenu); the
   resolved path used here matches the atlas registration. The
   menu-path mismatch is carried in frontmatter
   `unresolved_ambiguities` for operator confirmation.

4. In the MSA dialog, set the **Cluster** input to the column produced
   in step 2 (the `RandBetween(0, 5)` formula column).

### Scenario 2 — Alignment parameters button toggles dialog inputs

5. Click the **Alignment parameters** button in the MSA dialog. Verify
   the button adds input parameters to the dialog (the dialog grows
   the per-engine PepSeA alignment parameters surface — engine
   methods include `mafft --auto`, `mafft`, `linsi`, `ginsi`, `einsi`,
   `fftns`, `fftnsi`, `nwns`, `nwnsi` per atlas
   `bio.engines.msa-pepsea`). Click again to confirm the button
   toggles the parameter surface; the dialog returns to its prior
   shape.

### Scenario 3 — Run MSA via PepSeA and verify result column

6. Click **OK** to run alignment. MSA dispatches by alphabet:
   non-canonical (HELM) routes to the PepSeA Docker engine (atlas
   `bio.engines.msa-pepsea` — `package.ts#L1003`,
   `pepseaMsa` / role `sequenceMSA`). The Datagrok Docker subsystem
   manages PepSeA container lifecycle on demand.

7. Verify the result aligned column is added to the table (a new MSA
   column appears alongside the source HELM sequence column).

8. Verify the MSA column-header renderer is set on the result column
   (atlas `bio.rendering.column-header` — WebLogo header +
   conservation tracks paint on the aligned column; the renderer
   dispatches on the column's MSA flag).

9. Verify the per-cluster alignment invariant holds: for each Cluster
   value in the `RandBetween(0, 5)` column, the aligned sequences in
   that cluster are of the same length (PepSeA aligns within each
   Cluster group; cross-cluster lengths may differ).

## Notes

- **Source-text fixes silently applied** during migration (recorded in
  frontmatter `source_text_fixes`):
  - Step 5 had a trailing run of spaces after "properly"; trimmed.
  - Original step 6 said only "Check the new column." — vague. The
    assertion tail line ("Everything is good if the new MSA column is
    added, the renderer for the column is set, and the sequences
    within a cluster are of the same length.") has been promoted into
    three explicit verification steps (Scenario 3, steps 7-9 above)
    so the assertion surface is decidable. Per D-STEP-02 the
    expected-result line is preserved verbatim in semantic content;
    no assertion was dropped.
- **Sub-features covered:**
  - `bio.analyze.msa` (atlas L422) — feature root.
  - `bio.analyze.msa.dialog` (atlas L431) — `Bio | Analyze | MSA...`
    top-menu surface (`multipleSequenceAlignmentDialog`).
  - `bio.engines.msa-pepsea` (atlas L740) — non-canonical PepSeA
    Docker engine path (HELM is non-canonical → PepSeA via Docker
    container).
  - `bio.rendering.column-header` (atlas — MSA-aware column-header
    renderer) — exercised in step 8 verification.
  - `bio.detector` (atlas L122) — touched in setup (Macromolecule
    classification of the HELM sequence column on dataset open).
  Maps directly onto atlas critical path `bio.cp.msa-pepsea` (p1,
  `derived_from: public/packages/Bio/src/package.ts#L1003`).
- **Related bugs** surfaced for downstream awareness; per-bug
  cross-cutting coverage is delegated to chain-level
  `bug_focused_candidates`:
  - GROK-15176 — Bio's to-atomic-level produces molfiles with invalid
    isotope on heavy atoms (affects `bio.analyze.msa` +
    `bio.analyze.msa.dialog` as upstream producers). Chain proposes
    `bio-grok-15176-spec.ts` spanning `msa.md:Step 3`, `pepsea.md:Step
    3`, `convert.md:Step 2`; the cross-cutting MSA → To Atomic Level
    invariant is not asserted here. See chain
    `bug_focused_candidates[GROK-15176]`.
- **Cross-package contract** (atlas
  `bio.x.docker-container-eviction-msa-fallback`, coverage_type
  `edge`): PepSeA Docker container is on-demand; if the container is
  evicted between MSA invocations, the next non-canonical MSA must
  restart the container (await `DockerContainerStatus` transition)
  without surfacing an opaque crash. This scenario exercises the
  happy-path PepSeA dispatch; container-eviction recovery is
  addressed by the proactive lifecycle spec
  `bio-lifecycle-pepsea-container-spec.ts`
  (chain `proactive_lifecycle_specs[pepsea_container]`,
  `external_deps: [DockerContainer]`).
- **Helper coverage:** `bio.flow.msa` already registered in
  `helpers-registry.yaml:219` →
  `.claude/skills/grok-browser/references/bio.md:139`. The single
  registered helper covers both `msa.md` (canonical counterpart) and
  `pepsea.md` (this scenario). No new candidate helpers proposed.
- **Counterpart scenarios:** `msa.md` covers the canonical
  (FASTA → kalign WASM) branch on the same MSA top-menu surface;
  atlas `bio.cp.msa-canonical` (p0) is its canonical critical path.
  Together `msa.md` + `pepsea.md` realize both atlas-declared MSA
  critical paths.
- **Unresolved ambiguities** (carried in frontmatter
  `unresolved_ambiguities`):
  - `menu-path-bio-msa-vs-bio-analyze-msa` — chain analyzer
    (`bio.yaml` ambiguity for `pepsea.md:Step 3`) noted the source
    text opens MSA directly under "Bio" (no Analyze submenu); atlas
    registers MSA under "Bio | Analyze | MSA..." (`package.ts#L968`).
    Either the scenario predates the Analyze-submenu nesting or it
    was authored against a different menu version. This scenario uses
    the atlas-aligned `Bio > Analyze > MSA...` path; operator
    confirmation of the current menu path is needed.
  - `cluster-column-input-type-contract-integer-vs-categorical` —
    source does not specify whether the MSA dialog's Cluster input
    accepts an integer column (such as the `RandBetween(0, 5)`
    column produced in step 2) directly, or whether it expects a
    categorical/string column. This scenario follows the source text
    and supplies the integer column as-is; operator clarification of
    the Cluster-input column-type contract is needed before Automator
    can encode a typed assertion. (Same ambiguity also flagged on
    `msa.md`.)
- **`derived_from:` provenance:** the atlas entries cited above were
  derived from code anchors only — `bio.analyze.msa.dialog` and
  `bio.cp.msa-pepsea` both derive from `package.ts#L968` /
  `package.ts#L1003`. No help-doc derivations were used (per the
  binding sourcing rule).

---
{
  "order": 7
}
