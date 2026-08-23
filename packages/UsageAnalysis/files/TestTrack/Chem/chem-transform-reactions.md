---
feature: chem
realizes_atlas:
  - chem.cp.transform-reactions
realizes: []
priority: p1
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
related_bugs: []
realized_as:
  - chem-transform-reactions-spec.ts
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-transform-reactions-r4
    timestamp: 2026-08-22T13:46:00Z
    spec_runs:
      - spec: chem-transform-reactions-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 21
        failure_keys: []
        run_mode: headless-cold
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-11
    timestamp: "2026-08-20T19:30:00Z"
    failure_keys: []
    review_round: 1
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-transform-reactions-r3
    timestamp: 2026-08-22T13:52:14Z
    failure_keys: []
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      A new column named Desalted(canonical_smiles) (or equivalent auto-name) is
      added to the table with semType Molecule. Each value is the canonical
      SMILES of its input with water, hydrogen-halide, halide and sulfate
      fragments dropped: no output value still carries one of those fragments,
      no output value has more disconnected fragments than its input, and at
      least one output value differs from its input. Rows that were already a
      single fragment keep the same molecule — their canonical form is
      unchanged.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      The Desalted column cells render as molecule structures in the grid
      (semType Molecule is set and the grid binds the Molecule cell renderer to
      the column) and no cell is empty for rows whose input was a valid SMILES.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      A product column is added to the table. For the ester-hydrolysis reaction
      SMARTS [C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[O:3] applied to smiles.csv,
      at least one product cell holds a SMILES that differs from the desalted
      input on that row — the reaction ran and produced a carboxylic acid — and
      the number of such rows is strictly less than the number of valid input
      rows. Rows where the SMARTS finds no match carry the unchanged desalted
      input rather than an empty cell, so no valid-input row is left empty.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      The reaction-completion balloon reports a product count strictly greater
      than zero and strictly less than the total row count, and that count
      equals the number of product cells that differ from their desalted input —
      the reported number and the column agree.
  - anchor: "Scenario 3 Step 1"
    expectation: >-
      deprotect_fmoc.csv opens with its smiles column typed as Molecule, and a
      substructure search for the dialog's shipped default fragment
      O=C(N[*:1])OCC1c2ccccc2-c2ccccc21 (attachment point dropped, because RDKit
      matches a dummy atom only against another dummy atom) finds it on exactly
      five of the ten rows — the five Fmoc-protected amino acids. A different
      count means the fixture no longer matches the shipped fragment.
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      A deprotected(smiles) column (or similarly auto-named) is added with
      semType Molecule and no row yields an empty cell. Comparing input and
      output canonically — a matched row comes back as a molblock, an unmatched
      one as the input SMILES — exactly the five matching rows change and the
      five rows without a protecting group pass through as the same molecule. No
      output still contains the protecting group.
  - anchor: "Scenario 4 Step 4"
    expectation: >-
      A two-component product column is added. The Amide Coupling reaction
      SMARTS [C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3] is applied pairwise
      to two genuinely different reactant sets — the smiles1 column holds
      carboxylic acids, the smiles2 column holds primary amines — so each
      product is the amide of two distinct molecules. The ten acid/amine rows
      each yield a product and the two negative-control rows (one with no acid,
      one with no amine) yield none, so the non-empty product count is exactly
      ten, and every non-empty product differs from both reactant values on its
      row.
  - anchor: "Scenario 4 Step 5"
    expectation: >-
      All non-empty cells in the product column render as molecule structures in
      the grid, confirming semType Molecule is set on the column and the grid
      binds the Molecule cell renderer to it.
---

# Chem — Transform: Reactions surface

## Setup

1. Open Datagrok in a browser and wait for the platform to finish loading.
2. Open Demo Files > chem > smiles.csv and wait until the table opens and the canonical_smiles
   column renders molecule cells (semType Molecule detected, cells show structural drawings
   rather than raw text). The automated run loads this table programmatically rather than
   through File | Open — file-open actuation is not what this scenario covers.

## Scenarios

### Scenario 1: Remove Water and Salts strips multi-fragment SMILES to their largest organic fragment

Steps:

1. With smiles.csv open, choose the top menu Chem | Transform | Reactions | Remove Water
   and Salts.
2. In the Remove Water and Salts dialog, confirm that the table input is set to the open
   smiles.csv table and the molecules column is set to canonical_smiles (the detected Molecule
   column). Click OK.
3. Wait for the operation to complete. Observe that a new column appears in the grid
   alongside the original canonical_smiles column.
4. (Scenario 1 Step 4) Verify the new column. Every value is the canonical form of its input
   with water, hydrogen-halide, halide and sulfate fragments removed: no output value still
   carries one of those fragments, no output value has more disconnected fragments than its
   input, and at least one output value differs from its input. Rows that were already a
   single fragment describe the same molecule as before — only the SMILES text may have been
   rewritten into canonical form.
5. (Scenario 1 Step 5) Confirm the new column is typed as Molecule and that the grid draws
   its cells with the molecule cell renderer rather than as text.

Expected:
- The Desalted column is added with non-empty cells for all rows whose input was valid.
- Rows carrying a water or salt fragment lose it; no row gains a fragment.
- Rows that were already a single fragment describe the same molecule as before.
- The Desalted column is bound to the molecule cell renderer.

### Scenario 2: Transformation reaction applies a one-reactant SMARTS and writes a product column

Steps:

1. With smiles.csv still open, choose the top menu Chem | Transform | Reactions |
   Transformation.
2. The Transformation dialog opens with a reaction browser panel on the left and column /
   options inputs on the right. Click "+ New Reaction" and enter the ester-hydrolysis SMARTS
   [C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[O:3] in the reaction editor, then save it and select
   the new reaction card in the browser. Confirm the molecules column is set to
   canonical_smiles. Leave the Remove Salts option at its default (checked).
3. Click OK. Wait for the operation to complete.
4. (Scenario 2 Step 4) Verify that a product column appears in the grid and inspect its cells.
   At least one row must hold a product SMILES that differs from the desalted input on that
   row — the ester was hydrolysed to the free acid — and fewer rows than the total must differ,
   because not every molecule carries an ester. Rows the SMARTS did not match carry the
   unchanged desalted input, not an empty cell.
5. (Scenario 2 Step 5) Read the reaction-completion balloon. The product count it reports must
   be strictly greater than zero, strictly less than the total row count, and equal to the
   number of product cells that differ from their desalted input.

Expected:
- A product column is added with a non-empty cell on every valid input row.
- At least one product cell differs from its desalted input, and not all of them do.
- Unmatched rows carry the desalted input unchanged.
- The balloon's product count agrees with the number of changed product cells.

### Scenario 3: Deprotect removes the selected protecting group from matching rows

Steps:

1. (Scenario 3 Step 1) Open Files > App Data > Chem > tests > deprotect_fmoc.csv. This fixture
   holds one Molecule column, smiles, with ten rows: rows 1-5 are Fmoc-protected amino acids
   (Fmoc-Gly, Fmoc-Ala, Fmoc-Val, Fmoc-Leu, Fmoc-Phe) and rows 6-10 are the same acids plus
   benzoic acid with no protecting group. Before running anything, search the column for the
   dialog's default fragment and confirm it is found on exactly five rows — smiles.csv and the
   other bundled chemical fixtures contain no Fmoc carbamate at all, so on them every assertion
   below would be satisfied by the feature doing nothing.
2. With that table active, choose the top menu Chem | Transform | Reactions | Deprotect. The
   Deprotect dialog opens with a Molecules column selector and a Fragment field pre-filled with
   the default fluorenylmethyl-carbamate (Fmoc) protecting group
   O=C(N[*:1])OCC1c2ccccc2-c2ccccc21. Confirm the Molecules column is set to smiles and leave
   the Fragment field at its default value. Click OK.
3. (Scenario 3 Step 3) Wait for the operation to complete. A new column (named
   deprotected(smiles) or similar) is added to the grid, typed as Molecule, with no empty cell.
4. Verify the new column by comparing each output with its input as molecules rather than as
   text — a row the fragment matched comes back as a molblock, a row it did not comes back as
   the input SMILES, so a plain string comparison against a SMILES fragment can never see the
   product. Exactly the five matching rows must change, the five rows without a protecting group
   must describe the same molecule as before, and no output may still contain the protecting
   group.

Expected:
- The default fragment is found on exactly the five protected rows of the fixture.
- The deprotected column is added and typed as Molecule.
- Exactly those five rows change; the fragment is gone from every output.
- The five rows without the fragment pass through as the same molecule.

### Scenario 4: Two-Component Reaction combines reactants from two columns into a product column

Steps:

1. Open Files > App Data > Chem > tests > amide_coupling_2_columns.csv. This fixture carries two
   Molecule columns holding two genuinely different reactant sets: smiles1 holds ten simple
   carboxylic acids (benzoic, acetic, propionic, butyric, cyclohexanecarboxylic, p-toluic,
   4-chlorobenzoic, phenylacetic, nicotinic, 2-furoic) and smiles2 holds ten simple primary
   amines (aniline, benzylamine, cyclohexylamine, butylamine, p-toluidine, 4-chloroaniline,
   phenethylamine, 3-aminopyridine, ethylamine, p-anisidine). No acid also carries a primary
   amine and no amine also carries a carboxylic acid, so a product can only come from coupling
   the two columns. Two further rows are negative controls: row 11 holds toluene in place of an
   acid and row 12 holds anisole in place of an amine, so neither reacts.
2. With that table active, choose the top menu Chem | Transform | Reactions |
   Two-Component Reaction.
3. The Two-Component Reaction dialog opens with a reaction browser panel, two column selectors
   (Reactant 1 and Reactant 2), a Combination Mode selector and a Remove Salts toggle. Confirm
   Reactant 1 is smiles1, Reactant 2 is smiles2 and Combination Mode is pairwise, then select
   the pre-defined Amide Coupling card in the reaction browser
   ([C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]). Click OK.
4. (Scenario 4 Step 4) Wait for the operation to complete. A product column appears. Each of the
   ten acid/amine rows yields the amide of its two reactants; the two negative-control rows
   yield nothing. At least one product cell must be non-empty, fewer than all rows may carry a
   product, and every non-empty product must differ from both reactant values on its row —
   with two different reactant sets this is a real claim about a coupled product rather than a
   comparison against one molecule twice.
5. (Scenario 4 Step 5) Confirm the product column is typed as Molecule and that the grid draws
   its cells with the molecule cell renderer.

Expected:
- The Two-Component Reaction dialog opens with both reactant selectors and the reaction browser.
- After running, the acid/amine rows carry products and the two negative-control rows do not.
- Non-empty product cells contain SMILES strings that differ from both source reactants.
- The product column is bound to the molecule cell renderer.

## Automation notes

- The table is loaded through the file API rather than through File | Open; file-open
  actuation is covered elsewhere in this section.
- Scenario 3 reads its fixture as System:AppData/Chem/tests/deprotect_fmoc.csv.
- Scenario 3's five-matching-rows check guards the fixture, not the product's shipped default
  fragment: the substructure query is derived from the fragment constant the spec mirrors from
  Chem/src/analysis/deprotect.ts, so a changed shipped default leaves the count at five and this
  check silent. What fails on a changed default is the changed-rows comparison in Scenario 3
  step 3-4 — a dialog defaulting to another fragment leaves the Fmoc rows untouched, so the set of
  changed rows no longer equals the five matched rows. The query drops the [*:1] attachment point,
  because RDKit matches a dummy atom only against another dummy atom and the fragment as written
  matches nothing.
- Open follow-up: the Deprotect dialog's Fragment field is never read, so nothing observes the
  shipped default directly — the scenario only infers it from which rows changed. Reading that
  field would need a selector for it, and no browser reference doc names one; adding it requires
  live UI reconnaissance.
- Scenario 4 reads its fixture as System:AppData/Chem/tests/amide_coupling_2_columns.csv.
- Scenario 4 reads the Reactant 1, Reactant 2 and Combination Mode inputs rather than setting
  them: the fixture's two Molecule columns are what the dialog picks by default, and pairwise is
  the default mode. Selecting a reaction card and clicking OK are real UI gestures.
- Scenario 2's changed-cell count is compared against the Desalted column produced in
  Scenario 1, because the transformation desalts and canonicalizes each reactant before
  running the SMARTS — comparing against the raw input column would count desalting as a
  reaction product.
