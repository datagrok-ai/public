---
feature: chem
realizes_atlas:
  - chem.cp.calculate-clustering
realizes: []
priority: p1
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
ui_companion:
  - chem-calculate-clustering-ui.md
realized_as:
  - chem-calculate-clustering-spec.ts
related_bugs: []
scope_reductions:
  - id: SR-01
    check: E-VACUOUS-ASSERT-01
    rationale: |
      Scenario 3 does not exercise a computed maximum common substructure, and
      no assertion in it claims to.
      The Cluster input of the Cluster MCS dialog is registered
      `type: 'categorical'` (public/packages/Chem/src/package.ts#L833), which
      resolves to `Column.isCategorical`
      (core/shared/ddt/lib/src/data_frame/column_filter.dart#L56) — true only for
      string and bool columns. The BitBIRCH cluster-id column produced by
      Scenario 1 is an int column
      (public/packages/Chem/src/analysis/bit-birch/bitbirch-clustering.ts#L37,
      #L58) and is therefore never offered in that selector. In
      `System:DemoFiles/chem/smiles.csv` the only categorical column is
      `canonical_smiles` itself, and its 1000 values are 1000 distinct strings.
      So the partition the dialog pre-selects puts every molecule in a cluster
      of one, and the multi-member branch of `clusterMCS`
      (public/packages/Chem/src/rdkit-service/rdkit-service.ts#L633-L648) never
      runs. No categorical column with repeats is reachable in this scenario,
      by fixture or by any earlier step, so the multi-member path is left
      uncovered rather than faked.
      What is uncovered: that the value written for a cluster of two or more
      molecules is a substructure common to those molecules — the substance of
      the feature.
      What S3.4 asserts instead: the column is appended, is tagged
      `semType: Molecule` (set only on the success path,
      `package.ts#L860`), the row count is unchanged, and every row carries a
      non-empty structure string — which fails on the empty-entry failure paths
      (pre-filled '' at `rdkit-service.ts#L618`, the 45s worker-timeout restart
      at #L634-L645, the per-worker catch at #L643-L645, and a failed
      `get_smiles()` re-parse at #L667).
      An earlier revision asserted "at least one cell differs from that row's
      own molecule" and claimed such a cell could only come from a real
      multi-member computation. That is false and must not be restored: every
      entry, singleton pass-through included, is re-emitted through
      `getQueryMolSafe(it).get_smiles()`
      (`rdkit-service.ts#L659-L669`), so a passed-through molecule returns
      RDKit-canonicalized and differs from the source string whenever the
      source spelling was not already canonical. The 2026-08-22 run measured
      exactly that shape: 1000 rows, 929 differing, 71 already canonical.
      The premise itself is guarded: S3.4 asserts the largest cluster is
      exactly one member, so if a fixture change or a change to the categorical
      filter ever makes a multi-member partition reachable, this reduction
      fails loudly instead of going stale in silence — a promise about what is
      not covered has to break when it stops being true.
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-calculate-clustering-r4
    timestamp: "2026-08-22T15:41:00Z"
    spec_runs:
      - spec: chem-calculate-clustering-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 18
        failure_keys: []
        run_mode: headless-cold
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-08
    timestamp: "2026-08-19T10:00:00Z"
    failure_keys: []
    review_round: 2
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-calculate-clustering-r2
    timestamp: 2026-08-22T15:47:10Z
    failure_keys: []
expected_results:
  - anchor: "S1.4"
    expectation: >-
      A BitBIRCH cluster-id column is appended whose distinct-value count is
      greater than 1 and strictly less than the row count; the row count is
      unchanged.
  - anchor: "S3.4"
    expectation: >-
      A Cluster MCS column is appended, tagged as a Molecule column, with every
      row carrying a non-empty structure string; the row count is unchanged.
      Nothing here establishes that a common substructure was computed — see
      scope reduction SR-01.
  - anchor: "S4.4"
    expectation: >-
      A similarity-matrix table is produced with one similarity column per row
      (one per molecule), whose diagonal cells all read 1.0 and at least one
      off-diagonal cell differs from 1.0.
---

# Chem — Calculate Clustering (BitBIRCH, Butina, Cluster MCS, Similarity Matrix)

## Setup

1. Open Datagrok in a clean browser session and sign in.
2. Open `System:DemoFiles/chem/smiles.csv` via Browse > Demo Files > Chem > smiles.csv.
   Wait until the Grid viewer is visible and the Molecule column cells render structural
   images (not raw SMILES strings). Record the row count shown in the table status bar
   as the baseline row count.
   Actuation channel: a manual run opens the file through Browse; the paired spec opens it
   programmatically (`grok.dapi.files.readCsv`) because file-open actuation is not the
   subject of this scenario — it is covered by `chem-import-export-formats`. Everything
   from the Chem top menu onward is driven through the UI in both channels.

## Scenarios

### Scenario 1: BitBIRCH Clustering — cluster-id column with non-trivial partition

Steps:

1. **S1.1** — With `smiles.csv` open and molecule cells rendered, navigate to the top menu
   and select Chem > Calculate > BitBIRCH Clustering. The BitBIRCH Clustering dialog opens,
   showing a molecule-column selector (auto-detected) and at minimum a Threshold control.
2. **S1.2** — Confirm the molecule column selector shows the detected Molecule column.
   Accept all default values (threshold, fingerprint type) without changing them.
3. **S1.3** — Click OK. A progress indicator may appear briefly in the taskbar while the
   WASM BitBIRCH engine processes the molecules.
4. **S1.4** — Observe the grid. A cluster-id column (its name contains "Cluster" or "cluster") is
   appended as the last column. Verify:
   - The number of distinct values in the cluster-id column is greater than 1 (multiple
     clusters formed, not all molecules placed in one cluster).
   - The number of distinct values is strictly less than the row count (not every molecule
     assigned its own unique cluster).
   - The baseline row count is unchanged.
   - No JavaScript console errors fired during dialog open, OK click, or column-append.

### Scenario 2: Butina Cluster — delegated to the manual companion

Butina Cluster and the BitBIRCH ↔ Butina contrast are **not run here**. They live in
`chem-calculate-clustering-ui.md` (`target_layer: manual-only`), which carries the full
steps and expectations.

Reason for the delegation: Butina is a server-side Python script. Live recon on dev
(2026-08-19) opened the dialog and dispatched OK cleanly, but the job never joined its
`cluster (Butina)` column within ~90 s even on a 50-row extract — no error, no balloon, no
column. Asserting the column in an automated run would produce a false RED (the trigger
fires but the observable never lands), so the scenario is executed by hand instead.

Run the companion when you execute this file manually.

### Scenario 3: Cluster MCS — the column is produced, on the singleton partition only

Steps:

1. **S3.1** — From the top menu, select Chem > Calculate > Cluster MCS. The Cluster MCS
   dialog opens, showing a **Molecules** column selector and a required **Cluster** column
   selector (the categorical column that assigns each molecule to a cluster).
2. **S3.2** — In the Molecules selector, confirm the Molecule column is selected. Leave the
   Cluster selector on its default. Do not change either selector.
   The Cluster selector offers only categorical (string or bool) columns, so the cluster-id
   column produced in Scenario 1 — an int column — is not among them. In `smiles.csv` the only
   categorical column is `canonical_smiles` itself, and its 1000 values are all distinct, so
   the partition this run uses puts every molecule in a cluster of one. Note the name of the
   column the dialog pre-selected and how many clusters it yields; both are recorded in the run.
3. **S3.3** — Click OK and wait for the operation to complete.
4. **S3.4** — A Cluster MCS column (its name contains "MCS", "mcs", or "Scaffold") is appended to the
   grid. Verify:
   - Every row of this column holds a non-empty SMILES string. Empty cells are the product's
     failure signature: entries start empty, and stay empty when the worker times out, errors,
     or fails to re-parse the result.
   - The column is tagged as a Molecule column (its cells render as structures, not as raw
     text). The tag is applied only once the substructure computation returns, so an untagged
     column is a run that failed and was swallowed.
   - The baseline row count is unchanged.
   - No JavaScript console errors fired.

   None of this shows that a common substructure was *computed*: on a singleton partition each
   cell is that row's own molecule passed through and re-canonicalized. Scope reduction SR-01 in
   the frontmatter states what stays uncovered and why no reachable partition would cover it.

### Scenario 4: Similarity Matrix — diagonal reads as self-similarity

Precondition: the similarity matrix is pairwise, so a full `smiles.csv` (1000 rows) would
build a 1000x1000 matrix. Before running this scenario, take the first 50 molecules of
`smiles.csv` into a separate table named `clustering_matrix_subset` and make it the active
table; the matrix is computed over that bounded subset. (The paired spec builds this table
with `DG.DataFrame.create` / `DG.Column.fromStrings` and re-types the column as Molecule.)

Steps:

1. **S4.1** — With `clustering_matrix_subset` active, select Chem > Calculate > Similarity
   Matrix from the top menu. The Similarity Matrix dialog opens, showing a molecule-column
   selector.
2. **S4.2** — Confirm the molecule column is auto-detected. Accept all default values
   (fingerprint type, similarity metric). Click OK.
3. **S4.3** — Wait for the operation to complete; the pairwise run may take several seconds.
4. **S4.4** — A similarity-matrix result table appears among the open tables (its name
   contains "Similarity Matrix"): one row per molecule plus one similarity column per
   molecule. Verify:
   - The number of similarity columns equals the number of rows — one column per molecule.
     Without that, column *i* is not molecule *i* and the "diagonal" below is not the diagonal.
   - The diagonal values (the similarity of each molecule to itself) are all 1.0.
   - At least one off-diagonal cell differs from 1.0, confirming that pairwise comparisons
     between distinct molecules were computed and are non-trivial.
   - No JavaScript console errors fired during the run or during viewer rendering.

## Automation notes

- target_layer rationale: Scenarios 1, 3 and 4 require multi-step UI interaction with Chem
  top-menu dialogs, and the key assertions (column presence, distinct-value count, diagonal
  values in a grid viewer) are DOM- and data-state assertions that Playwright drives natively.
  Scenario 2 (Butina) is manual-only and lives in `chem-calculate-clustering-ui.md`.
- Distinctness assertion: read the cluster-id column's category count and assert it is greater
  than 1 and less than the total row count. Do not assert an exact count — the threshold default
  may produce different partition counts on different engine versions.
- Column-appearance waits: snapshot the column list **before** clicking OK, not after. A
  snapshot taken after the trigger can already contain the new column, in which case the
  wait can never be satisfied, degenerates into a sleep for its full timeout, and — with the
  timeout swallowed — still passes on the follow-up probe. That makes run duration unstable
  and is a flake source, not merely lost time. Assert on the column name the wait returns.
- Console errors: the trap is installed at test start (before the first dialog opens) and
  named ambient noise is excluded by name, never by a broad `Failed to load resource` pattern
  that could hide a real failure.
- What a failed Cluster MCS run actually looks like: there is no error cell. The per-cluster
  result array is pre-filled with empty strings
  (`public/packages/Chem/src/rdkit-service/rdkit-service.ts#L618`); a cluster of one molecule
  gets that molecule copied through (`#L632`) and then canonicalized like every other entry
  (`#L659-L669`); the 45-second worker timeout restarts
  the worker and leaves the entry empty (`#L634-L645`); and the one failure handler raises the
  balloon `Cluster MCS Error` plus a console line and returns the column untouched
  (`public/packages/Chem/src/package.ts#L864-L866`, column created at `#L850`). A run that fails
  outright therefore produces an all-empty, untagged column — never a cell holding an error
  message. An earlier version of this scenario asked the tester to confirm the absence of a raw
  error cell: that is a state the product has no code path to produce, so the check could not
  fail and proved nothing. The observables that do distinguish success from failure are the
  Molecule tag on the column and the absence of empty cells; both are what S3.4 now asks for.
  Do not reinstate the error-cell check.
  A cell differing from its own row's molecule is **not** such an observable, and an earlier
  revision of this file wrongly said it was. Every entry — singleton pass-through included — is
  re-emitted through `getQueryMolSafe(it).get_smiles()` (`rdkit-service.ts#L659-L669`), so a
  passed-through molecule comes back RDKit-canonicalized and differs from the source string
  whenever the source spelling was not already canonical. Do not restore that check either.
- derived_from: atlas entry `chem.cp.calculate-clustering` references
  `public/packages/Chem/src/package.ts#L797` (BitBIRCH registration) and
  `public/help/datagrok/solutions/domains/chem/scripts/butina-cluster.md` for the Butina script.
