---
feature: bio
sub_features_covered:
  - bio.transform.to-atomic-level
  - bio.transform.to-atomic-level.action
  - bio.transform.to-atomic-level.single
  - bio.transform.to-atomic-level.api
  - bio.transform.helm-to-mol
  - bio.transform.molecules-to-helm
target_layer: playwright
coverage_type: regression
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

Gate F SR extension scenario (cycle 2026-06-01-bio-migrate-02) that
realizes atlas critical path `bio.cp.to-atomic-level` (p1) and
cross-feature interaction `bio.x.bio-to-chem-via-atomic-level`
(`coverage_type: regression`, `related_bugs: [GROK-15176]`). Targets
five net-new bio sub-features over the live covered union and one
already-covered umbrella id (`bio.transform.to-atomic-level`) carried
for density and as the dispatch anchor for the cell-action /
top-menu paths.

The scenario exercises the four entry points into the bio→molfile→helm
conversion surface: (a) `Bio | Transform | To Atomic Level...` top-menu
on a Macromolecule column, (b) the cell action `to atomic level`
surfaced from the column-header context-menu (atlas
`bio.transform.to-atomic-level.action`), (c) the public API wrappers
`seq2atomic(seq, nonlinear)` (atlas `bio.transform.to-atomic-level.api`)
and `toAtomicLevelSingleSeq(sequence)` (atlas
`bio.transform.to-atomic-level.single`), and (d) the round-trip path
`Bio | Transform | Molecules to HELM...` consuming the molfile output
back into HELM (atlas `bio.transform.molecules-to-helm`). The
`getMolFromHelm` helper-converter pipeline (atlas
`bio.transform.helm-to-mol`, `package.ts#L1684`) is the underlying
service exercised by the To-Atomic-Level path on HELM input — verified
via the V3K molfile column shape that the cell renderer dispatch
consumes downstream.

## Setup

Datasets (canonical Bio test fixtures, two notation classes):

- `System.AppData/Bio/tests/filter_HELM.csv` — HELM notation. Required
  for the `helm-to-mol` pipeline (atomic conversion routes through the
  HELM-to-molfile converter under `src/utils/helm-to-molfile/converter/`).
- `System.AppData/Bio/tests/filter_FASTA.csv` — FASTA notation. Routes
  through the linear `_toAtomicLevel` branch (atlas
  `bio.transform.to-atomic-level` description, `package.ts#L863`).

Entry points (atlas labels):

- **Bio | Transform | To Atomic Level...** — top-menu action,
  `toAtomicLevel` function (`package.ts#L863`).
- **Bio | Transform | Molecules to HELM...** — top-menu action,
  `moleculesToHelmTopMenu` (`package.ts#L824`).
- **right-click sequence column → to atomic level** — cell action
  `toAtomicLevelAction` (`package.ts#L886`); atlas
  `bio.transform.to-atomic-level.action`.
- `await grok.functions.call('Bio:seq2atomic', { seq, nonlinear })` —
  public API wrapper (`package.ts#L1605`); atlas
  `bio.transform.to-atomic-level.api`.
- `await grok.functions.call('Bio:toAtomicLevelSingleSeq', { sequence })` —
  single-sequence variant (`package.ts#L911`); atlas
  `bio.transform.to-atomic-level.single`.

## Scenarios

### Scenario 1 — Top-menu To Atomic Level on HELM input

1. Open `System.AppData/Bio/tests/filter_HELM.csv` from the Files
   browser. The Macromolecule detector classifies the sequence column
   synchronously (atlas `bio.detector`); the column is tagged with
   `quality=Macromolecule, units=helm`.
2. On the menu ribbon open **Bio | Transform | To Atomic Level...**
   The To-Atomic-Level dialog appears with `name=dialog-To-Atomic-Level`
   and the HELM column prefilled in the column input.
3. Accept the prefilled column and click **OK**. The pipeline routes
   through `HelmToMolfileConverter` (atlas `bio.transform.helm-to-mol`,
   `package.ts#L1684`).

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
2. Click **to atomic level** in the menu (atlas
   `bio.transform.to-atomic-level.action`, registered as a cell action;
   `package.ts#L886`).

Expected:
- The same molfile column shape from Scenario 1 appears; the cell
  action runs the same `toAtomicLevel` code path against the clicked
  Macromolecule column (no separate dialog when the column is already
  bound by the action target).
- No console error; the new column is appended to the active
  dataframe and renders via the Chem cell renderer (downstream
  cross-package contract from interaction
  `bio.x.bio-to-chem-via-atomic-level`).

### Scenario 3 — Top-menu To Atomic Level on FASTA input (linear branch)

1. Open `System.AppData/Bio/tests/filter_FASTA.csv` (linear FASTA
   sequences; detector tags `units=fasta`).
2. On the menu ribbon open **Bio | Transform | To Atomic Level...**
   The dialog appears with the FASTA Macromolecule column prefilled.
3. Click **OK** to run the linear `_toAtomicLevel` branch (atlas
   `bio.transform.to-atomic-level` description, `package.ts#L863`).

Expected:
- A new `molfile(Sequence)` column appears with `semType=Molecule`,
  `units=molblock`. The branch difference vs Scenario 1 is internal
  (linear vs HELM-pipeline routing); the column-level contract is the
  same. This step seeds the GROK-15176 regression guard: the produced
  molfile MUST NOT contain isotope flags on heavy atoms (see Scenario
  5 cross-check).

### Scenario 4 — Molecules to HELM round-trip

1. With the molfile column from Scenario 1 (or Scenario 3) present on
   the active dataframe, open **Bio | Transform | Molecules to HELM...**
   (atlas `bio.transform.molecules-to-helm`, `package.ts#L824`). The
   dialog appears with the molfile column prefilled.
2. Click **OK** to run the Python-script + monomer-library-matching
   converter that produces a HELM column from molecule input.

Expected:
- A new column appears whose `semType` is `Macromolecule` and whose
  cell values are HELM strings (e.g. `PEPTIDE1{...}$$$$V2.0`).
- The round-trip yields a HELM column structurally consumable by the
  Bio detector when the column is re-detected — verifying the molfile
  output of To-Atomic-Level is a valid input to the
  `moleculesToHelm` reverse pipeline (the
  `bio.x.bio-to-chem-via-atomic-level` cross-feature contract is
  exercised end-to-end).

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

- target_layer rationale: the top-menu and cell-action surfaces
  (`bio.transform.to-atomic-level`, `.action`) are UI-only triggers
  whose dispatch into `toAtomicLevel` runs only via a real menu /
  context-menu pointer event; an `apitest`-only scenario cannot
  exercise the action-registration path that GROK-15176's repro
  context describes (cell-action route from MSA-column header).
  Playwright is the dominant layer; Scenario 5 carries the API
  wrappers and the regression guard inside the Playwright runner
  as `grok.functions.call(...)` calls (no separate apitest needed).
- Atlas citations:
  - `bio.cp.to-atomic-level` (p1, `bio.yaml#L1349`) — open HELM
    dataset → Bio | Transform | To Atomic Level... → V3K molfile
    column produced → Chem renderer renders cells. Scenarios 1+2
    realize this critical path on HELM; Scenario 3 covers the linear
    FASTA branch.
  - `bio.x.bio-to-chem-via-atomic-level` (`bio.yaml#L1476`,
    `coverage_type: regression`, `related_bugs: [GROK-15176]`).
    Scenario 4's molecules-to-HELM round-trip and Scenario 5's
    isotope-flag regression guard cover the cross-package
    output-validity contract this interaction names.
  - `edge_cases[].source_bug: GROK-15176` (`bio.yaml#L1699`) —
    `sub_features: {bio.analyze.msa, bio.analyze.msa.dialog,
    bio.transform.to-atomic-level, .api, .action}`. Scenario 5's
    explicit `MASS=1` heavy-atom assertion against both wrapper
    outputs is the targeted regression guard for this edge case.
- Net-new ids vs live_covered_union (51 ids on incoming chain
  rev 14 + bio-renderer-dispatch.md = 51): five of the six ids on
  this scenario are net-new — `.action`, `.single`, `.api`,
  `helm-to-mol`, `molecules-to-helm`. The umbrella id
  `bio.transform.to-atomic-level` is already in the union via
  `convert.md`'s sub_features; carried here for density and to
  identify the dispatch surface the .action and top-menu paths
  bind to. Projected non-manual_only union after this scenario:
  49 + 5 = 54 of 99 (~54.5%); still below the 70% threshold —
  remaining first-pass and second-pass proposals in the chain
  `gaps[].proposed_action` block continue the breadth lift.
- Bug-aware tie-break (STEP C): selected over equivalently-sized
  first-pass candidates because GROK-15176 is one of the chain's
  five `bug_focused_candidates[]` entries (the
  `bio-grok-15176-spec.ts` proposal that does not yet exist on
  disk). Per the tightened F-BUG-COVERAGE-01 branch (i), this
  scenario's `related_bugs: [GROK-15176]` is the same closure
  evidence `convert.md`, `msa.md`, and `pepsea.md` already
  provide; the chain's downstream `bug-grok-15176-spec.ts`
  proposal remains a valid Automator seed for the cross-cutting
  Chem-renderer / PubChem integration test, orthogonal to this
  scenario.
- Deferrals (none): every sub_feature on this scenario is reached
  by at least one of the five scenarios above; no
  technical-dependency blocker requires deferral.

---
{
  "order": 21
}
