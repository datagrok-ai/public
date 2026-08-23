---
feature: chem
realizes_atlas:
  - chem.cp.calculate-properties
realizes: []
priority: p0
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
realized_as:
  - chem-calculate-properties-spec.ts
related_bugs: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-08
    timestamp: "2026-08-19T10:00:00Z"
    failure_keys: []
    review_round: 2
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-calculate-properties-r3
    timestamp: 2026-08-22T13:56:05Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-calculate-properties-r2
    timestamp: 2026-08-22T13:32:10Z
    spec_runs:
      - spec: chem-calculate-properties-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 89
        run_mode: headless-cold
        failure_keys: []
expected_results:
  - anchor: "S1.4"
    expectation: >-
      An MW column (or the selected property column) is appended to the active
      table grid. The column contains numeric values greater than zero for valid
      SMILES rows. The total grid row count is unchanged. No JavaScript console
      errors fire during dialog open, OK click, or column-append.
  - anchor: "S1.6"
    expectation: >-
      Nine property columns are appended to the grid. Eight carry the fixed
      names HBA, HBD, LogP, LogS, PSA, Rotatable bonds, Stereo centers and
      Molecule charge; the ninth is MW re-appended under a de-duplicated name
      because an MW column already exists from S1.4. Every one of the nine
      carries a numeric value on every row — smiles.csv holds no unparseable
      molecule, which S3.3 establishes independently by producing an InChI
      string on all rows. The grid row count is unchanged. No JavaScript console
      errors fire during the dialog re-open, OK click, or column-append.
  - anchor: "S2.4"
    expectation: >-
      A Mutagenicity column is appended to the grid. It is a string column whose
      values come from the OCL risk-level vocabulary — "Unknown", "None", "Low",
      "High" — and it is non-empty on every molecule row. At least one row must
      read "Low" or "High": smiles.csv holds drug-like molecules, so a column
      that carries no assessed risk means the risk engine produced no result.
      "Unknown" is the level OCL writes when it failed to assess the molecule,
      so it is not accepted as a result — nor is "None".
  - anchor: "S2.6"
    expectation: >-
      All four toxicity-risk columns (Mutagenicity, Tumorigenicity, Irritating
      effects, Reproductive effects) are appended. Each is a string column from
      the same "Unknown"/"None"/"Low"/"High" vocabulary, non-empty on every row,
      and each carries at least one row reading "Low" or "High". "Unknown" is
      the failed-to-assess level and does not count as a result. The grid row
      count is unchanged.
  - anchor: "S3.3"
    expectation: >-
      A column named `inchi` is appended to the grid. Every row carries a string
      starting with "InChI=" (e.g. "InChI=1S/..."), and no row that held a
      molecule has a blank or null value.
  - anchor: "S4.3"
    expectation: >-
      A column named `inchi_key` is appended to the grid. Every row carries a
      27-character InChIKey matching [A-Z]{14}-[A-Z]{10}-[A-Z] (e.g.
      "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"), no molecule row is blank, and every value
      differs from the `inchi` value on the same row.
---

# Chem — Calculate Properties, Toxicity Risks, To InChI, To InChI Keys

## Setup

1. Open Datagrok in a clean browser session and sign in.
2. Open `System:DemoFiles/chem/smiles.csv` via Browse > Demo Files > Chem > smiles.csv.
   Wait until the Grid viewer is visible and its canvas is attached. Record the table's row
   count as the baseline row count, and confirm the `canonical_smiles` column is typed as
   Molecule (semType `Molecule`) — that is the readiness signal the paired spec waits on.
   A manual executor reads the row count from the table status bar and confirms the Molecule
   cells render as structures rather than raw SMILES text.
   Actuation channel: a manual run opens the file through Browse; the paired spec opens it
   programmatically (`grok.dapi.files.readCsv`) because file-open actuation is not the
   subject of this scenario — it is covered by `chem-import-export-formats`. Everything
   from the Chem top menu onward is driven through the UI in both channels.

## Scenarios

### Scenario 1: Chemical Properties — MW only, then full property set

Steps:

1. **S1.1** — With `smiles.csv` open and molecule cells rendered, navigate to the top menu
   and select Chem > Calculate > Chemical Properties. The Chemical Properties dialog opens,
   showing a **Molecules** column selector (auto-detected) and a list of property checkboxes
   (MW, HBA, HBD, LogP, LogS, PSA, Rotatable bonds, Stereo centers, Molecule charge).
2. **S1.2** — In the dialog, confirm the Molecules selector shows `canonical_smiles`.
   Ensure only the MW checkbox is selected (check MW, uncheck all others). This produces
   the minimal single-property case.
3. **S1.3** — Click OK.
4. **S1.4** — Observe the grid. An `MW` column of type double is appended as the last
   column, holding values greater than zero. The baseline row count from Setup is unchanged.
5. **S1.5** — Navigate to Chem > Calculate > Chemical Properties again. The dialog re-opens
   with the same molecule column selected.
6. **S1.6** — In the dialog, select all nine property checkboxes: MW, HBA, HBD, LogP, LogS,
   PSA, Rotatable bonds, Stereo centers, Molecule charge. Click OK. Nine new columns are
   appended to the right of the existing ones. Eight of them carry the fixed names HBA, HBD,
   LogP, LogS, PSA, Rotatable bonds, Stereo centers and Molecule charge; the ninth is MW
   calculated again and appended under a **de-duplicated name** (`MW (2)` or similar),
   because the `MW` name is already taken by the column from S1.4. That is why the
   fixed-name checklist has eight entries and not nine — it is not a missing column. Every
   one of the nine carries a numeric value on every row. The grid row count equals the
   baseline.

### Scenario 2: Toxicity Risks — Mutagenicity only, then full risk set

Steps:

1. **S2.1** — From the top menu, select Chem > Calculate > Toxicity Risks. The Toxicity
   Risks dialog opens with the **Molecules** selector on `canonical_smiles` and four risk
   checkboxes: Mutagenicity, Tumorigenicity, Irritating effects, Reproductive effects.
2. **S2.2** — Confirm only Mutagenicity is checked (it is the default initial value per the
   source registration). Leave the other three unchecked. Click OK.
3. **S2.3** — The dialog closes.
4. **S2.4** — A `Mutagenicity` column is appended to the grid. It is a **string** column
   holding OpenChemLib risk levels — one of `Unknown`, `None`, `Low`, `High` — never a
   number or a boolean. It is non-empty on every row, and at least one row reads `Low` or
   `High`: smiles.csv holds drug-like molecules, so a column with no assessed risk means the
   risk engine produced no result. `Unknown` is the level OCL writes when it failed to assess
   the molecule; it is not accepted as a result, and neither is `None`. The vocabulary is
   closed at those four values with no intermediate level, so requiring `Low` or `High` cannot
   fail on a healthy product. The grid row count equals the baseline.
5. **S2.5** — Navigate to Chem > Calculate > Toxicity Risks again.
6. **S2.6** — In the dialog, select all four checkboxes: Mutagenicity, Tumorigenicity,
   Irritating effects, Reproductive effects. Click OK. Four risk columns are appended, each
   from the same `Unknown`/`None`/`Low`/`High` vocabulary, each non-empty on every row and
   each with at least one row reading `Low` or `High` — `Unknown` is the failed-to-assess
   level and does not count as a result. Mutagenicity is recalculated and lands
   under a de-duplicated name, as in S1.6, because the `Mutagenicity` name is taken by the
   column from S2.4. The grid row count equals the baseline.

### Scenario 3: To InchI — column of InChI strings appended

Steps:

1. **S3.1** — From the top menu, select Chem > Calculate > To InchI... (that is the menu
   label as registered). A function-call dialog opens.
2. **S3.2** — Click OK. Wait for the operation to complete (a progress indicator may appear
   briefly in the taskbar).
3. **S3.3** — A column named `inchi` is appended to the grid. Every row shows an InChI
   string starting with "InChI="; no molecule row is blank or null. The grid row count
   equals the baseline.

### Scenario 4: To InchI Keys — column of InChI Key strings appended

Steps:

1. **S4.1** — From the top menu, select Chem > Calculate > To InchI Keys.... A function-call
   dialog opens.
2. **S4.2** — Click OK and wait for the operation to complete.
3. **S4.3** — A column named `inchi_key` is appended to the grid. Every row shows a
   27-character InChIKey matching [A-Z]{14}-[A-Z]{10}-[A-Z] (for example,
   BSYNRYMUTXBXSQ-UHFFFAOYSA-N); no molecule row is blank, and every value differs from the
   `inchi` value on the same row. The grid row count equals the baseline.

## Automation notes

- Grid column presence: after each OK click, query the table's column list (e.g.
  `df.columns.names()` via JS API in the browser context or by asserting the column
  header DOM element is visible) to confirm the expected column name was appended.
- Row count guard: assert `df.rowCount === baselineRowCount` after each calculator
  completes. A row count change indicates the calculator wrote a new table instead of
  appending columns.
- For Scenarios 3 and 4 the top-menu action always opens a function-call dialog
  (`[name="dialog-To-InchI"]` / `[name="dialog-To-InchI-Keys"]`). Wait for its OK button,
  click it, then wait for the appended column. There is no immediate-execution branch.
- Actuation channel: the table is opened through the JS API (`grok.dapi.files.readCsv`),
  as declared in Setup. From there on nothing is injected — the whole menu path
  (top-menu click > dialog > checkbox > OK) is driven through the UI, so the dialog's
  column detection and property selection are genuinely exercised. Do not substitute a
  `grok.functions.call` for the menu path.
- Column de-duplication: the second Chemical Properties run and the second Toxicity Risks
  run recalculate a property whose column name is already taken, so the platform appends it
  under a de-duplicated name. Assert on the fixed names that are new (eight, and three) plus
  the total number of appended columns (nine, and four) — never on nine or four fixed names.
