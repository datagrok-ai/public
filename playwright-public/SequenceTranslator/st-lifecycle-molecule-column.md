# SequenceTranslator — Molecule Column Lifecycle

## Setup

1. Ensure the SequenceTranslator package is loaded.
2. Open Datagrok.
3. Load a table containing a Molecule (SMILES) column. A fixture with
   Markush core structures (e.g., `chem_enum_cores.csv` or any table with
   a SMILES column containing R-group labels such as `[*:1]`) can be used.
   If not available, create a DataFrame with a column named `smiles` containing
   at least two SMILES strings with R-group placeholders.

## Scenarios

### Scenario 1: Enumerate Markush Structure from top-menu on Molecule column

Steps:
1. Select a Molecule cell in the table (a SMILES string with R-group
   labels such as `[*:1]`, `[*:2]`).
2. Navigate to `Chem | Transform | Markush Enumeration...`.
3. In the `polyToolEnumerateChemUI` dialog, the core SMILES should be
   seeded from the selected cell.
4. Configure one R-group list (e.g., `[CH3, C2H5, C3H7]`) and choose
   `Cartesian` mode.
5. Click Enumerate.

Expected:
- The dialog opens seeded with the selected Molecule cell's SMILES.
- The enumeration runs in Cartesian mode.
- A result molecule table is produced with row count equal to the
  Cartesian product of cores × R-group lists.
- All result rows are valid SMILES strings (canonicalized via a batched
  `Chem:convertNotation` call).

### Scenario 2: Enumerate Markush Structure from cell context menu

Steps:
1. Right-click a Molecule grid cell (semType = `Molecule`) in the table.
2. Select `Enumerate Markush Structure...` from the context menu.
3. The `polyToolEnumerateChemUI` dialog opens seeded with the clicked cell.
4. Accept defaults and click Enumerate.

Expected:
- The context menu item `Enumerate Markush Structure...` is present on
  right-click of a Molecule cell (wired via `detectors.js` context-menu
  dispatch for Molecule semType).
- The dialog is seeded with the clicked cell's SMILES value.
- The enumeration completes and a result column is produced.

### Scenario 3: Zip mode vs Cartesian mode enumeration

Steps:
1. Select a Molecule cell with two R-group positions (`[*:1]` and `[*:2]`).
2. Open `Chem | Transform | Markush Enumeration...`.
3. Configure two R-group lists of equal length.
4. Choose `Zip` mode; click Enumerate.
5. Note the row count of the result.
6. Reopen the dialog for the same cell; switch to `Cartesian` mode.
7. Click Enumerate; note the row count.

Expected:
- Zip mode produces row count equal to the length of one R-group list
  (pairwise substitution).
- Cartesian mode produces row count equal to
  `len(list1) × len(list2)`.
- All result rows are valid SMILES.