---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [bio.cp.msa-pepsea]
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

Covers the non-canonical Multiple Sequence Alignment path: runs the
PepSeA Docker engine on a HELM (non-canonical peptide alphabet)
dataset with a per-row Cluster column, then verifies the aligned
result column is added, the MSA column-header renderer is wired, and
sequences within the same cluster come out the same length.

The canonical (FASTA → kalign WASM) counterpart lives in `msa.md`;
this scenario covers the PepSeA Docker-engine branch on HELM.

## Setup

Dataset: `System.AppData/Bio/tests/filter_HELM.csv` (HELM macromolecule
fixture for this section). The dataset is opened so the Macromolecule
detector classifies the sequence column synchronously — the
`Bio | Analyze | MSA...` top menu is only reachable when a
Macromolecule column is present on the active table view, and
HELM-aligned data routes MSA through the non-canonical PepSeA Docker
engine instead of the canonical kalign WASM path.

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
   dialog opens (`multipleSequenceAlignmentDialog`). Source text
   listed the path as **Bio** > **MSA** (no Analyze submenu); the
   path used here is the current one.

4. In the MSA dialog, set the **Cluster** input to the column produced
   in step 2 (the `RandBetween(0, 5)` formula column).

### Scenario 2 — Alignment parameters button toggles dialog inputs

5. Click the **Alignment parameters** button in the MSA dialog. Verify
   the button adds input parameters to the dialog (the dialog grows
   the per-engine PepSeA alignment parameters surface — engine
   methods include `mafft --auto`, `mafft`, `linsi`, `ginsi`, `einsi`,
   `fftns`, `fftnsi`, `nwns`, `nwnsi`). Click again to confirm the
   button toggles the parameter surface; the dialog returns to its
   prior shape.

### Scenario 3 — Run MSA via PepSeA and verify result column

6. Click **OK** to run alignment. MSA dispatches by alphabet:
   non-canonical (HELM) routes to the PepSeA Docker engine
   (`pepseaMsa` / role `sequenceMSA`). The Datagrok Docker subsystem
   manages PepSeA container lifecycle on demand.

7. Verify the result aligned column is added to the table (a new MSA
   column appears alongside the source HELM sequence column).

8. Verify the MSA column-header renderer is set on the result column
   (WebLogo header + conservation tracks paint on the aligned
   column; the renderer dispatches on the column's MSA flag).

9. Verify the per-cluster alignment invariant holds: for each Cluster
   value in the `RandBetween(0, 5)` column, the aligned sequences in
   that cluster are of the same length (PepSeA aligns within each
   Cluster group; cross-cluster lengths may differ).

## Notes

- This scenario exercises the happy-path PepSeA dispatch only. If
  the PepSeA container gets evicted between MSA invocations, the
  next run must restart it automatically rather than crash — that
  recovery path is covered separately by
  `bio-lifecycle-pepsea-container.md`.
- The molfile isotope-flag validity issue for GROK-15176 is
  exercised more thoroughly in `bio-transform-atomic-level.md`.
- Open question: it's not clear whether the MSA dialog's Cluster
  input is meant to accept a raw integer column (like the
  `RandBetween(0, 5)` column used here) or expects a
  categorical/string column instead — this scenario follows the
  original test text and supplies the integer column as-is (same
  open question as `msa.md`).

---
{
  "order": 7
}
