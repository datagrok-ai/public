---
feature: chem
realizes_atlas:
  - chem.cp.calculate-identifiers
realizes: []
priority: p1
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
produced_from: atlas-driven
realized_as:
  - chem-calculate-identifiers-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers:
  - openChemMenuItem
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-08
    timestamp: "2026-08-20T06:00:00Z"
    failure_keys: []
    review_round: 2
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-calculate-identifiers
    timestamp: 2026-08-22T12:38:55Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-calculate-identifiers
    timestamp: 2026-08-22T12:14:40Z
    spec_runs:
      - spec: chem-calculate-identifiers-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 36
        failure_keys: []
        run_mode: headless-cold
expected_results:
  - anchor: "S1.5"
    expectation: >-
      An identifier column is appended to the active table grid. It is named
      after the chosen toSource value verbatim — `chembl`, not a display label
      such as "ChEMBL ID". The resolved-row count in the appended column is
      greater than zero: at least one input molecule maps to a known identifier
      in the chosen source, and an all-blank column means the resolution did not
      run. The baseline row count is unchanged.
  - anchor: "S1.6"
    expectation: >-
      A second identifier column from a different toSource is appended alongside
      the first, named `pubchem`. It too resolves at least one row. The values
      in the two columns differ from each other in at least one row, confirming
      that the two resolution passes ran against distinct external sources
      rather than copying the same result. The baseline row count is still
      unchanged.
  - anchor: "S2.4"
    expectation: >-
      The Chemical Properties calculator appends `MW` — the one property whose
      parameter is enabled by default — and appends nothing that is not one of
      its own nine properties (MW, HBA, HBD, LogP, LogS, PSA, Rotatable bonds,
      Stereo centers, Molecule charge). Each row of the appended column carries
      a numeric value — not an error string, not a raw stack trace, and not the
      original SMILES string passed through unchanged. The baseline row count is
      unchanged.
  - anchor: "S2.5"
    expectation: >-
      A second pass with a different calculator appends that calculator's own
      column (`clogP` for logP) and appends nothing else, so the two passes
      produced different property kinds rather than one calculator running
      twice. Every row of the new column holds a finite number. The baseline row
      count is unchanged.
  - anchor: "S3.4"
    expectation: >-
      A new table labelled "conformers" (or containing "conformer") opens in the
      workspace. Its "smiles" column reads `CCCCC` on every row — the molecule
      the step sketched, not the script's own default of `CCCC`; a table built
      from the default means the sketcher gesture never reached the run and
      leaves every other check below unable to fail. The table contains at least
      two rows (pentane with default settings produces multiple conformers). The
      "conformer" column holds incrementing integer values starting at 1. The
      "molblock" column holds non-empty molblock strings containing atom
      coordinates and "M  END". The "energy" column holds MMFF94 energies: this
      run leaves Optimize at its default of true, so every row must be a finite
      number — an all-NaN energy column means optimization silently failed. NaN
      would be acceptable only on a run with Optimize set to false, which this
      scenario does not perform. The "rmsd" column exists; the first row reads
      0.0 (reference conformer); at least one other row shows a non-zero RMSD.
      An all-empty conformers table or a single-row result is a failure, not a
      pass.
---

# Chem — Map Identifiers, Biochemical Properties, Generate Conformers

## Setup

1. Open Datagrok in a clean browser session and sign in.
2. Open `System:DemoFiles/chem/smiles.csv` via Browse > Demo Files > Chem > smiles.csv.
   Wait until the Grid viewer is visible and the `canonical_smiles` cells render structural
   images (not raw SMILES strings). Record the row count shown in the table status bar
   as the baseline row count before any scenario runs.
   Actuation channel: a manual run opens the file through Browse; the paired spec opens it
   programmatically (`grok.dapi.files.readCsv`) because file-open actuation is not the
   subject of this scenario — it is covered by `chem-import-export-formats`. Everything
   from the Chem top menu onward is driven through the UI in both channels.

## Scenarios

### Scenario 1: Map Identifiers — resolved identifier column with non-zero hit count

Steps:

1. **S1.1** — With `smiles.csv` open and molecule cells rendered, open the top menu and
   select Chem > Calculate > Map Identifiers.... The Map Identifiers dialog opens. It shows
   a molecule-column selector (auto-detected to `canonical_smiles`), a "From" source
   selector, and a "To" source selector.
2. **S1.2** — In the dialog, confirm the molecule column selector shows the detected
   `canonical_smiles` column. In the "From" source selector, choose `smiles`. In the "To"
   source selector, choose `chembl`. The source entries are lowercase identifier-source ids
   (`smiles`, `chembl`, `pubchem`, `inchi`, ...), not display labels such as "ChEMBL ID".
3. **S1.3** — Click OK. A progress indicator appears in the taskbar while the identifier
   resolution runs against the chosen external source.
4. **S1.4** — Wait for the operation to complete (the progress indicator closes).
5. **S1.5** — Observe the grid. A column named `chembl` — the toSource value verbatim — is
   appended as the last column. Verify that it holds at least one non-blank resolved
   identifier value; an entirely blank or null column means the resolution did not run,
   which is a failure. The baseline row count is unchanged.
6. **S1.6** — Open Chem > Calculate > Map Identifiers... again. In the "To" source
   selector, choose `pubchem`. Click OK and wait for completion. A second column named
   `pubchem` is appended, itself resolving at least one row. Compare the two columns row by
   row: at least one row must show different values between them, confirming that the two
   resolution passes targeted distinct external databases.

### Scenario 2: Biochemical Properties — property columns from two independent calculators

Steps:

1. **S2.1** — From the top menu, select Chem > Calculate > Biochemical Properties (the menu
   label carries no "..." suffix). The Biochemical Properties dialog opens. It shows a table
   selector, a molecule-column selector, a navigation list of available biochemical
   calculators on the left (Chemical Properties, logD, logP, pI, pKa), and an editor panel
   for the currently selected calculator on the right.
2. **S2.2** — Confirm the molecule column selector shows the `canonical_smiles` column from
   `smiles.csv`. In the navigation list on the left, tick the checkbox on the **Chemical
   Properties** row. Confirm the editor panel shows the calculator's parameters with the
   molecule column referenced.
3. **S2.3** — Click OK. Wait for the calculator to complete (a progress indicator may appear
   briefly in the taskbar) and for the Biochemical Properties dialog to close.
4. **S2.4** — Observe the grid. Chemical Properties appends `MW` — the only one of its nine
   properties whose parameter is enabled by default — and appends nothing outside that set
   of nine (MW, HBA, HBD, LogP, LogS, PSA, Rotatable bonds, Stereo centers, Molecule
   charge). Each row shows a numeric value: not an error string, not a raw stack trace, and
   not the original SMILES string passed through unchanged. The baseline row count is
   unchanged.
5. **S2.5** — Open Chem > Calculate > Biochemical Properties again. In the navigation list,
   tick **logP** instead. Click OK and wait for completion. That calculator appends its own
   column, named `clogP`, and appends nothing else — a different property kind from the
   first pass's `MW`, which is what confirms two independent calculators ran rather than one
   running twice. Every row of `clogP` holds a finite number. The baseline row count is
   unchanged.

### Scenario 3: Generate Conformers — conformer table with 3D molblock output

Steps:

1. **S3.1** — From the top menu, select Chem > Calculate > Generate Conformers.... The
   Generate Conformers script dialog opens. It shows a single-molecule input (a sketcher
   canvas — this script takes one molecule at a time), Num conformers, Optimize, RMS
   threshold, Max attempts and Random seed. Read the six inputs and confirm they stand at
   50, true (checked), 0.1, 5000 and 42; the rest of the scenario depends on those defaults
   holding, and a changed Optimize default would invalidate the finite-energy expectation in
   S3.4.
2. **S3.2** — Click the molecule input to open the sketcher, enter the SMILES string `CCCCC`
   (pentane — a flexible molecule with rotatable bonds that produces multiple distinct
   conformers) and confirm the sketcher. Read the sketcher's molecule field back before
   confirming: it must show `CCCCC`. The molecule field of the script dialog is pre-populated
   with the script's own default of `CCCC` (butane), so entering butane would leave the run
   indistinguishable from one where the sketcher gesture never landed; the whole scenario
   would still pass on the default. Leave all other parameters at their default values.
3. **S3.3** — Click OK. A progress indicator appears in the taskbar while the server-side
   Python script (RDKit ETKDGv3) runs.
4. **S3.4** — The script completes and a new table labelled "conformers" (or containing
   "conformer") opens in the workspace as a new tab. Verify:
   - The "smiles" column holds the molecule that was sketched — `CCCCC` on every row, not the
     script's default `CCCC`. A table built from the default means the sketcher gesture never
     reached the run, and every other check below would pass regardless.
   - The table has at least two rows. A single-row result indicates conformer generation
     found only one conformer, which is a failure for a flexible molecule like pentane.
   - The "conformer" column contains incrementing integer values starting at 1.
   - The "molblock" column contains non-empty molblock strings — each cell contains atom
     coordinates and ends with "M  END", not an empty string or the original SMILES.
   - The "energy" column contains MMFF94 energies. Optimize was left at its default of
     true, so every row must be a finite number; an all-NaN energy column means the
     optimization silently failed. NaN would be acceptable only on a run with Optimize set
     to false, which this scenario does not perform.
   - The "rmsd" column exists. The first row shows 0.0 (the reference conformer has zero
     RMSD relative to itself). At least one other row shows a non-zero RMSD, confirming
     the generated conformers differ geometrically from the reference.
   - An all-empty conformers table or a zero-row result is a failure, not a pass.

## Automation notes

- target_layer rationale: all three scenarios require multi-step UI dialogs driven via
  the Chem top menu, and the assertions (column presence, non-empty resolved values,
  non-zero hit count for identifiers; conformer table row count and molblock content)
  are DOM- and data-state assertions that Playwright drives natively. Generate Conformers
  is a server-side Python script, but its observable output is a new Datagrok table with
  accessible column values — no special Playwright support is needed beyond polling for
  the new table to appear.
- coverage_type rationale: this scenario walks the happy path of three calculators to prove
  they still produce their documented output. It does not reproduce a reported defect, so it
  is `smoke` rather than `regression`, and `related_bugs` is empty by construction.
- openChemMenuItem helper: use the `openChemMenuItem` helper (registered
  `helpers/chem.ts:3`) to navigate nested Chem top-menu leaves. The Chem menu can be
  hidden by viewport overflow; the helper dispatches a click event that survives that.
- Map Identifiers: the from/to source selectors are rendered inside a custom editor
  widget (`MapIdentifiersEditor`) as plain `<select>` elements addressed by
  `[name="input-From-Source"]` / `[name="input-To-Source"]` inside
  `[name="dialog-Map-Identifiers"]`. Drive them with a real `selectOption` gesture — a bare
  `select.value = …` assignment leaves the `ui.input.choice` model at its default and the
  run appends a wrongly-named column. After OK, poll for a new column named after the
  chosen toSource with a timeout of 90 seconds.
- Biochemical Properties: the navigator list items are `.biochem-calc-nav-item` rows inside
  the `biochem-calc-nav-list` container, each an `input[type=checkbox]` plus a span label.
  Tick the checkbox of the calculator the scenario names — never "the first visible item",
  which does not pin which calculator ran. A real Playwright `check()` reaches the Dart
  model; a synthetic `change` event does not. Poll for the new column with a timeout of
  120 seconds.
- Generate Conformers: the molecule input is a `.d4-input-molecule-canvas`, not a text
  field; clicking it opens a sketcher dialog carrying an `input[placeholder*="SMILES"]`.
  Fill that input, press Enter, then click the sketcher's OK — all as real gestures. After
  the outer OK, poll for a new table whose name contains "conformer" in `grok.shell.tables`.
  Access its "conformer", "molblock", "energy", and "rmsd" columns to run the structured
  assertions. Timeout: 180 seconds (server-side Python script with MMFF94 optimization).
- Do not pre-populate inputs via JS API — the full menu path must be driven through the
  UI so the dialog's column-detection and parameter-selection logic is exercised. Every
  control this scenario touches was confirmed on dev (2026-08-22) to respond to a real
  Playwright gesture, so there is no exception to this rule in this pair.
- derived_from: atlas entry `chem.cp.calculate-identifiers` references
  `public/packages/Chem/src/package.ts#L691` (Map Identifiers registration) and
  `public/packages/Chem/scripts/generate_conformers.py#L6` (Generate Conformers
  top-menu registration).
