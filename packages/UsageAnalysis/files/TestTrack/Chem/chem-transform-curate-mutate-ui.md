---
feature: chem
realizes_atlas: []   # untyped: no matching atlas type
realizes:
  - chem.cp.transform-curate-mutate
priority: p2
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  Curate (Chem | Transform | Curate...) and Mutate (Chem | Transform | Mutate...)
  are server-side Python scripts that dispatch to the Chem Python compute service.
  Live MCP recon on dev.datagrok.ai 2026-08-20 (cycle 2026-08-19-chem-new-09):
  the Curate dialog opened with the exact defaults the scenario names and OK
  dispatched the job (status bar "Running Curate..."), but curated_molecule never
  joined the table within ~6 minutes — no column, no error, no balloon. A direct
  Chem:Mutate call on the scenario's default SMILES did not return within 110s
  either. Every expected-results bullet depends on that compute output, so a
  Playwright spec would time out at service dispatch (a false-RED, not a product
  defect). Executed by hand: run each scenario, wait for the compute to land, and
  read the curated_molecule / mutations values per the Expected: clauses. See
  migration-lessons.md section 2 (backend/service-dependent features -> manual-only)
  and the sibling chem-calculate-clustering-ui.md Butina disposition.
related_bugs:
  - id: GROK-20747
    status: open
---

# Chem — Transform: Curate and Mutate

## Setup

1. Open the demo dataset DemoFiles > chem > chem_standards.csv. Verify the
   table loads with 14 rows, a Name column, and a smiles column whose cells
   render as molecule structures (semType Molecule detected).

## Scenarios

### Scenario 1: Curate normalises targeted rows and leaves others unchanged

Steps:
1. With chem_standards.csv open, click the top menu Chem > Transform > Curate.
2. In the Curate dialog, confirm the molecules column is pre-selected in the
   "molecules" field. Leave all boolean options at their defaults (normalization
   ON, reionization ON, neutralization ON, mainFragment ON, kekulization OFF,
   tautomerization OFF). Click OK.
3. Wait for the script to complete. Observe that a new column named
   curated_molecule appears appended to the table.
4. Verify the curated_molecule column: for the row named metal_non
   (input SMILES CCC(=O)O[Na]) the curated value must differ from the input
   (the metal-containing salt is stripped to its organic fragment); for the row
   named metal_st (input CCC(=O)[O-].[Na+]) the curated value must also differ.
   For rows already in the canonical standardised form (e.g. parent_st, norm_st)
   the curated value must equal the input SMILES exactly — confirming the script
   is selective, not a blanket transformer.
5. Re-run Curate with only normalization=ON and all other options OFF. Confirm
   that the new curated_molecule column differs from the input only on the rows
   norm_non (input CN(C)=CC=C[O-], expected output includes the normalised
   zwitterion form) and reion_non, while all other rows are unchanged.

Expected:
- A curated_molecule column is added after the script completes.
- Rows targeted by the active option carry SMILES values that differ from the
  input; rows not targeted carry values identical to the input.
- An output column where every cell equals its input row is a test failure — it
  means the normalisation option did not apply at all.

### Scenario 2: Mutate generates parseable non-parent structures

Steps:
1. With any table open (chem_standards.csv is sufficient), click the top menu
   Chem > Transform > Mutate.
2. In the Mutate dialog, enter the SMILES CN1C(CC(O)C1=O)C1=CN=CC=C1 in the
   molecule field (the default from the script). Leave steps=1, randomize=true,
   maxRandomResults=100. Click OK.
3. Wait for the script to complete. A new table named mutations opens (or is
   added) with a mutations column. Verify: the row count is greater than zero
   and at most 100; every cell in the mutations column contains a non-empty
   SMILES string that parses as a valid molecule (the cell renders a structure,
   not an error icon); no cell in the mutations column is equal to the input
   SMILES CN1C(CC(O)C1=O)C1=CN=CC=C1.
4. Re-run Mutate with the same molecule and steps=2. Confirm: the row count is
   still greater than zero and at most 100; the mutations column contains no
   copy of the original parent SMILES; the structures in the table visibly
   differ from those returned in the steps=1 run (a two-step mutation space
   extends beyond the one-step space).

Expected:
- The mutations table is created with a non-empty mutations column.
- Every mutations cell parses as a valid molecule and differs from the parent.
- Increasing steps from 1 to 2 does not reduce the result to zero rows.

### Scenario 3: Curate kekulization option produces Kekule output

Steps:
1. Open chem_standards.csv. Click Chem > Transform > Curate.
2. In the Curate dialog, set kekulization=ON and leave all other options at
   their defaults (normalization ON, reionization ON, neutralization ON,
   mainFragment ON). Click OK.
3. After the curated_molecule column appears, inspect the cells for rows whose
   input SMILES contains aromatic notation (lower-case atom letters such as
   the benzene ring in parent_non: O=C([O-])c1ccccc1). Verify that the
   curated_molecule value for those rows is written in Kekule form — upper-case
   carbon atoms with explicit alternating single/double bond characters — and
   that no lower-case aromatic atom letter appears in the curated output for
   those rows.

Expected:
- The curated_molecule column contains Kekule SMILES for rows whose input has
  aromatic notation when kekulization=ON.
- At least one cell (e.g. parent_non row) must visibly differ in notation from
  the input: input uses 'c1ccccc1', output uses the Kekule equivalent.

## Notes

- target_layer rationale: both Curate and Mutate are invoked through top-menu
  dialogs (Chem > Transform > Curate and Chem > Transform > Mutate) and their
  outputs land as new columns or new tables visible in the grid. Playwright can
  drive the dialog inputs and assert the resulting column values via DOM
  inspection.

---
{
  "order": 7,
  "datasets": ["System:DemoFiles/chem/smiles.csv"]
}
