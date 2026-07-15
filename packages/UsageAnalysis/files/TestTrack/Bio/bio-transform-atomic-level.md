---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [bio.cp.to-atomic-level]
realizes: [bio.transform.to-atomic-level]
produced_from: atlas-driven
related_bugs:
  - GROK-15176
realized_as:
  - bio-transform-atomic-level-spec.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T05:10:00Z
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
    timestamp: 2026-06-02T04:35:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:50:00Z
    spec_runs:
      - spec: bio-transform-atomic-level-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 189
        failure_keys: []
---

# Bio | Transform — To Atomic Level + HELM/molfile round-trip

Covers the Bio→molfile→HELM conversion surface via its four entry
points: the **Bio | Transform | To Atomic Level...** top-menu on a
Macromolecule column, the **to atomic level** cell action on the
column-header context menu, the public API wrappers `seq2atomic` and
`toAtomicLevelSingleSeq`, and the round-trip path **Bio | Transform |
Molecules to HELM...** that converts a molfile column back into HELM.

Scenario 5 also guards against GROK-15176: a produced molfile must
never carry an illegal isotope flag on a heavy atom, or downstream
PubChem standardization silently rejects it.

## Setup

Datasets (canonical Bio test fixtures, two notation classes):

- `System.AppData/Bio/tests/filter_HELM.csv` — HELM notation. Required
  for the `helm-to-mol` pipeline (atomic conversion routes through the
  HELM-to-molfile converter under `src/utils/helm-to-molfile/converter/`).
- `System.AppData/Bio/tests/filter_FASTA.csv` — FASTA notation. Routes
  through the linear `_toAtomicLevel` branch.

Entry points:

- **Bio | Transform | To Atomic Level...** — top-menu action,
  `toAtomicLevel` function (`package.ts#L863`).
- **Bio | Transform | Molecules to HELM...** — top-menu action,
  `moleculesToHelmTopMenu` (`package.ts#L824`).
- **right-click sequence column → to atomic level** — cell action
  `toAtomicLevelAction` (`package.ts#L886`).
- `await grok.functions.call('Bio:seq2atomic', { seq, nonlinear })` —
  public API wrapper (`package.ts#L1605`).
- `await grok.functions.call('Bio:toAtomicLevelSingleSeq', { sequence })` —
  single-sequence variant (`package.ts#L911`).

## Scenarios

### Scenario 1 — Top-menu To Atomic Level on HELM input

1. Open `System.AppData/Bio/tests/filter_HELM.csv` from the Files
   browser. The Macromolecule detector classifies the sequence column
   synchronously; the column is tagged with
   `quality=Macromolecule, units=helm`.
2. On the menu ribbon open **Bio | Transform | To Atomic Level...**
   The To-Atomic-Level dialog appears with `name=dialog-To-Atomic-Level`
   and the HELM column prefilled in the column input.
3. Accept the prefilled column and click **OK**. The pipeline routes
   through `HelmToMolfileConverter`.

Expected:
- A new column appears whose `semType` is `Molecule` and whose units
  metadata is `molblock`.
- The column name follows the pattern `molfile(<source>)`; on a
  repeated run within the same view, the new column lands as
  `molfile(<source>) (2)` (matches the convert-notation matrix
  convention from `convert.md`).
- Cell content for a non-null source row is a V3K molfile string
  (`M  V30 BEGIN CTAB` header present in the cell value when inspected
  via `df.col('<name>').get(0)`).

### Scenario 2 — Column-header cell-action `to atomic level` on HELM

1. With `filter_HELM.csv` open, right-click the sequence column header
   to surface the column context menu.
2. Click **to atomic level** in the menu.

Expected:
- The same molfile column shape from Scenario 1 appears; the cell
  action runs the same `toAtomicLevel` code path against the clicked
  Macromolecule column (no separate dialog when the column is already
  bound by the action target).
- No console error; the new column is appended to the active
  dataframe and renders via the Chem cell renderer.

### Scenario 3 — Top-menu To Atomic Level on FASTA input (linear branch)

1. Open `System.AppData/Bio/tests/filter_FASTA.csv` (linear FASTA
   sequences; detector tags `units=fasta`).
2. On the menu ribbon open **Bio | Transform | To Atomic Level...**
   The dialog appears with the FASTA Macromolecule column prefilled.
3. Click **OK** to run the linear `_toAtomicLevel` branch.

Expected:
- A new `molfile(Sequence)` column appears with `semType=Molecule`,
  `units=molblock`. The branch difference vs Scenario 1 is internal
  (linear vs HELM-pipeline routing); the column-level contract is the
  same. This step seeds the GROK-15176 regression guard: the produced
  molfile MUST NOT contain isotope flags on heavy atoms (see Scenario
  5 cross-check).

### Scenario 4 — Molecules to HELM round-trip

1. With the molfile column from Scenario 1 (or Scenario 3) present on
   the active dataframe, open **Bio | Transform | Molecules to HELM...**.
   The dialog appears with the molfile column prefilled.
2. Click **OK** to run the Python-script + monomer-library-matching
   converter that produces a HELM column from molecule input.

Expected:
- A new column appears whose `semType` is `Macromolecule` and whose
  cell values are HELM strings (e.g. `PEPTIDE1{...}$$$$V2.0`).
- The round-trip yields a HELM column structurally consumable by the
  Bio detector when the column is re-detected — verifying the molfile
  output of To-Atomic-Level is a valid input to the
  `moleculesToHelm` reverse pipeline.

### Scenario 5 — Public API wrappers: seq2atomic + toAtomicLevelSingleSeq + GROK-15176 isotope-flag regression guard

This scenario calls the two public API wrappers directly via
`grok.functions.call` (no dataframe round-trip needed), and asserts
the molfile-validity invariant that GROK-15176 documents (heavy atoms
must not carry `isotope=1` flags — the downstream PubChem
standardization rejection symptom that the bug names).

Steps:

1. Run the linear single-sequence wrapper:

   ```ts
   const seq = 'ACDEFGHIK';
   const molLinear: string = await grok.functions.call(
     'Bio:toAtomicLevelSingleSeq', { sequence: seq });
   ```

2. Run the public API wrapper for the (potentially) non-linear path:

   ```ts
   const molNonLinear: string = await grok.functions.call(
     'Bio:seq2atomic', { seq: 'PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0',
                        nonlinear: true });
   ```

3. For each returned molfile string, inspect the V3000 atom block:

   ```ts
   function heavyAtomIsotopeFlags(mol: string): string[] {
     // V3000 atom lines look like: "M  V30 <idx> <symbol> x y z 0 ..."
     // Isotope is encoded via "MASS=<n>" key on the atom line.
     return mol.split(/\r?\n/)
       .filter((l) => /^M  V30 \d+ [A-Za-z]+ /.test(l))
       .filter((l) => / MASS=1\b/.test(l))    // isotope=1 marker
       .filter((l) => ! / [HD] /.test(l));    // exclude H / D (legitimate isotope use)
   }
   const offendersLinear = heavyAtomIsotopeFlags(molLinear);
   const offendersNonLinear = heavyAtomIsotopeFlags(molNonLinear);
   ```

Expected:
- `molLinear` and `molNonLinear` are non-empty strings whose first
  non-blank line matches `/V3000/` (V3K molfile shape).
- `offendersLinear.length === 0` AND `offendersNonLinear.length === 0`
  — no heavy atom carries an `isotope=1` MASS flag. This is the
  GROK-15176 fixed-behavior assertion. A non-zero count means the
  bug has regressed: downstream PubChem standardization on the
  produced molfile will reject heavy atoms with illegal isotopes,
  breaking the Chem panel identifier lookup the
  `bio.x.bio-to-chem-via-atomic-level` interaction depends on.
- The two wrapper code paths (`toAtomicLevelSingleSeq` →
  `_toAtomicLevel`; `seq2atomic` → linear or HELM-pipeline branch
  depending on `nonlinear`) both produce equally clean molfiles for
  the same monomer alphabet on identical sequence input.

## Notes

- No dedicated cross-cutting Chem-renderer / PubChem regression spec
  for GROK-15176 exists yet; `convert.md`, `msa.md`, and `pepsea.md`
  also touch this bug at a lighter level.

---
{
  "order": 21
}
