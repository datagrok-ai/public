---
feature: chem
realizes_atlas:
  - chem.cp.transform-notation-roundtrip
realizes: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-09
    timestamp: "2026-08-20T06:00:00Z"
    failure_keys: []
    review_round: 2
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-transform-notation-r3
    timestamp: 2026-08-22T15:23:00Z
    spec_runs:
      - spec: chem-transform-notation-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 28
        failure_keys: []
        run_mode: headless-cold
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-transform-notation-r3
    timestamp: 2026-08-22T15:31:31Z
    failure_keys: []
priority: p0
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
realized_as:
  - chem-transform-notation-spec.ts
related_bugs:
  - id: GROK-17964
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      Overwrite is off, so a new column (e.g. canonical_smiles_molblock) is
      added beside the original rather than replacing it. Its values are
      non-empty V2000 molblock strings for every valid input row, and none of
      them is a raw SMILES string.
  - anchor: "Scenario 1 Step 5b"
    expectation: >-
      A second Convert Notation run back to smiles produces a column whose
      canonical SMILES match the canonical SMILES of the original
      canonical_smiles column row for row over a 30-row sample — identity
      preserved, notation changed. At least one row's stored text visibly
      differs between the original column and the molblock column.
  - anchor: "Scenario 2 Step 3a"
    expectation: >-
      Join is on, so a new column is added beside the original. Its values are
      V2000 molblocks, and comparing them against the default-coordinate
      molblock column produced in Scenario 1 shows the coordinate block
      differing on at least one row while the atom and bond counts stay
      identical on every compared row — the coordinates were recalculated and
      the atom-bond graph was not touched.
  - anchor: "Scenario 2 Step 3b"
    expectation: >-
      A subsequent Convert Notation from the recalculated column back to
      canonical SMILES produces values that match the original canonical SMILES
      row for row over a 30-row sample, confirming from the SMILES side that the
      atom-bond graph survived the coordinate recalculation.
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      A new column is added with SMILES values for the names that ChEMBL could
      resolve; each such value is a parseable SMILES string. The resolved-row
      count is greater than zero and strictly less than the total row count,
      every row whose Name cell is empty yields an empty SMILES cell rather than
      an arbitrary value, and the row named ZZZZ NOT A COMPOUND ZZZZ — a name
      that is present but cannot be resolved — also yields an empty SMILES cell
      rather than a structure.
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      The resolved SMILES column has semType Molecule and the grid binds the
      Molecule cell renderer to it, so resolved rows draw as structures. Rows
      that produced no SMILES hold an empty string, which the molecule renderer
      leaves as an empty cell.
---

# Chem — Transform notation roundtrip

## Setup

1. Open Datagrok in a browser and wait for the platform to finish loading.
2. Open Demo Files > chem > smiles.csv and wait until the table opens and the
   canonical_smiles column renders molecule cells. The automated run loads this table
   programmatically rather than through File | Open — file-open actuation is not what this
   scenario covers.

## Scenarios

### Scenario 1: Convert Notation preserves molecular identity across a notation round-trip

Steps:

1. With smiles.csv open and the canonical_smiles column visible, choose top menu
   Chem | Transform | Convert Notation.
2. In the Convert Notation dialog, set the Molecules column to canonical_smiles, set Target
   Notation to molblock, leave Overwrite unchecked so a new column is added, and click OK.
3. Wait for the operation to complete and for the new column to appear in the grid.
4. (Scenario 1 Step 4) Verify that the new column (e.g. canonical_smiles_molblock) is present
   next to the original — Overwrite was off — and that its cells contain non-empty V2000
   molblock strings rather than raw SMILES text.
5. In the same table, choose top menu Chem | Transform | Convert Notation again.
   - (Scenario 1 Step 5a) Set the Molecules column to the newly added molblock column, set
     Target Notation to smiles, and click OK.
   - (Scenario 1 Step 5b) When the second new column appears, compare its values to the
     original canonical_smiles column: the canonical forms must match row for row, and the
     stored text in the intermediate molblock column must visibly differ from the SMILES text.

Expected:
- The molblock column (or equivalent auto-named column) is added beside the original and
  contains non-empty V2000 molblock strings for every valid input row.
- Converting that molblock column back to smiles produces canonical SMILES values that match
  the canonical SMILES of the original column row for row, confirming that molecular identity
  is preserved while only the notation changed.
- At least one row shows visibly different stored text between the original column and the
  molblock column (i.e. the conversion is not a no-op).

### Scenario 2: Recalculate Coordinates changes coordinates without altering the atom-bond graph

Steps:

1. With the same smiles.csv table active, choose top menu
   Chem | Transform | Recalculate Coordinates.
2. In the Recalculate Coordinates dialog, set the Molecules column to canonical_smiles, choose
   Method CoordGen, leave Join checked so a new column is added beside the original, and click
   OK. Wait for completion.
3. When the new column appears:
   - (Scenario 2 Step 3a) Compare it against the molblock column produced in Scenario 1, which
     carries the default coordinates for the same molecules. The atom coordinates must differ
     on at least one row, while the atom and bond counts must be identical on every compared
     row.
   - (Scenario 2 Step 3b) Run Convert Notation on the recalculated column back to canonical
     smiles (top menu Chem | Transform | Convert Notation, source = the recalculated column,
     target notation = smiles) and compare the result to the original canonical_smiles column
     row for row.

Expected:
- The coordinate-recalculated column is added beside the original and holds V2000 molblocks.
- Its coordinates differ from the default-coordinate molblocks of Scenario 1 on at least one
  row, while the atom and bond counts are unchanged — the geometry moved, the graph did not.
- Converting the recalculated column back to canonical SMILES produces values that match the
  original canonical SMILES row for row, confirming the same thing from the SMILES side.

### Scenario 3: Names To Smiles resolves known names and leaves the rest empty

Steps:

1. Open the file System:AppData/Chem/tests/names_to_smiles.csv, which holds a Name column with
   a handful of drug names, one blank row, and one row carrying the deliberately
   non-chemical name ZZZZ NOT A COMPOUND ZZZZ, which no chemical registry can match.
2. With the Name column visible, choose top menu Chem | Transform | Names To Smiles. In the
   dialog, set the Names column to Name and click OK. Wait for the operation to complete; it
   calls the ChEMBL name-resolution service and may take several seconds.
3. (Scenario 3 Step 3) When the new SMILES column appears, count the non-empty resolved rows.
   The resolved-row count must be greater than zero and strictly less than the total row
   count, the row whose Name is blank must carry an empty SMILES value, and the row named
   ZZZZ NOT A COMPOUND ZZZZ must carry an empty SMILES value too.
4. (Scenario 3 Step 4) Confirm the resolved SMILES column is typed as Molecule and that the
   grid draws its cells with the molecule cell renderer, so resolved rows appear as structures
   and rows with no SMILES stay empty.

Expected:
- The Names To Smiles operation adds a column with semType Molecule whose non-empty cells
  contain parseable SMILES strings and which the grid draws with the molecule cell renderer.
- The resolved-row count is greater than zero and strictly less than the total row count.
- A row with a blank name yields a blank SMILES value rather than an arbitrary one.
- A row whose name is present but unresolvable yields a blank SMILES value as well — the
  resolver rejects the name instead of inventing a structure for it.

## Automation notes

- The tables are loaded through the file API rather than through File | Open; file-open
  actuation is covered elsewhere in this section.
- Dialog inputs (Molecules, Target Notation, Method, Names) are set on the dialog's input
  model rather than by clicking the controls; clicking OK is a real UI gesture. The Overwrite
  and Join booleans are read back from the same model and asserted, so a flipped default is
  reported as such rather than surfacing later as a missing column.
- Canonical-SMILES comparisons call `grok.functions.call('Chem:canonicalize', {molecule})` on
  each sampled value of both columns before comparing. Direct string equality on the raw
  column values is insufficient because SMILES input may not be in canonical form.
- The canonical comparisons run over the first 30 valid rows rather than all 1000: each value
  costs a round-trip through the canonicalizer, and 30 rows keep the two comparisons inside
  the step budget while still covering a mix of ring systems, stereocentres and salts.
- Names To Smiles: the fixture carries two non-resolving rows — the blank one and the row
  named ZZZZ NOT A COMPOUND ZZZZ — so both halves of the claim are checked: a blank name
  yields a blank SMILES, and a name that is present but unresolvable yields a blank SMILES
  rather than a structure. The unresolvable row is identified by its Name value, not by a
  row index, so the check survives future fixture edits; the step also asserts that exactly
  one such row exists, so a fixture that lost it fails instead of passing vacuously. The
  name is deliberately non-chemical so no ChEMBL update can ever make it resolve.
