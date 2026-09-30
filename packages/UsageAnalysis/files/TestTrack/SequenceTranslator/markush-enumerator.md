---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [molecule_column]
realizes: [chem.menu.transform.markush-enumeration]
realized_as: []
related_bugs: [GROK-20223, GROK-20372]
---

# SequenceTranslator — Markush Enumerator

Checks the Markush Enumerator on a Molecule column with R-labelled cores. A core is
taken from a cell, R-groups are imported from a table, and Cartesian and Zip modes
produce the expected number of molecules. The result table name is honoured, and
the app view keeps its content when the user switches away and back.

## Setup

- Open `System:AppData/SequenceTranslator/tests/chem_enum_cores.csv` (column `Core`,
  4 cores with `[*:n]` labels) and
  `System:AppData/SequenceTranslator/tests/chem_enum_rgroups.csv` (columns `R1`,
  `R2`, `R3`; R1 and R2 have 4 distinct values each). All these columns are detected
  as Molecule.
- The core used below, `C1C2CN([*:2])CC2CN1[*:1]`, is row 2 of `chem_enum_cores` and
  has two R-positions, R1 and R2. So Cartesian mode gives 4 × 4 = 16 molecules and
  Zip mode gives 4.

## Scenarios

### Block A — Core from the cell context menu, R-groups from a table, Cartesian

1. In `chem_enum_cores`, right-click the cell `C1C2CN([*:2])CC2CN1[*:1]` and choose
   **Enumerate Markush Structure...**.
   * Expected result: the **Markush Enumerator** dialog opens with exactly one core
     card in **Cores**, subtitled `core 1 · R1, R2`. **Enumerator type** is
     **Cartesian**. **OK** is disabled.
2. In **R-Groups**, click the import (folder) icon. In **Import R-Groups** set
   **Table** = `chem_enum_rgroups`, **Column** = `R1`, **Target R#** = `1`, and click
   **OK**. Repeat with **Column** = `R2`, **Target R#** = `2`.
   * Expected result: **R-Groups** lists R1 and R2 with 4 substituents each. The
     status line next to **Enumerator type** reads **16 molecules will be
     generated**. **OK** is enabled.
3. Set **Table name** to `BddMarkushCartesian` and click **OK**.
   * Expected result: a table view **BddMarkushCartesian** opens with **16** rows and
     the columns **Enumerated**, **Core**, **R1**, **R2** (all rendered as molecules).
     No error balloon.

### Block B — Zip mode and table name (GROK-20223)

1. Reopen the dialog as in Block A step 1 and import R1 and R2 as in step 2.
2. Switch **Enumerator type** to **Zip**.
   * Expected result: the status line reads **4 molecules will be generated**.
3. Set **Table name** to `BddMarkushZip` and click **OK**.
   * Expected result: a table view **BddMarkushZip** opens with **4** rows.

### Block C — Top menu and app view (GROK-20372)

1. Make the `Core` cell of row 2 of `chem_enum_cores` current and open **Chem |
   Transform | Markush Enumeration...**.
   * Expected result: the **Markush Enumerator** dialog opens with that core in
     **Cores**. Close it with **CANCEL**.
2. Open the **Markush Enumerator** app (**Browse > Apps > Chem > Markush
   Enumerator**).
   * Expected result: a view **Markush Enumerator** shows the **Cores**,
     **R-Groups** and **Preview** areas; the ribbon has the **Enumerate** button and
     the **Enumerator type** selector.
3. Click **Enumerate** in the ribbon.
   * Expected result: an **Enumerate** dialog opens with **Output** = **New table**,
     **Table name**, **Remove duplicates** and a **Run** button. Close it with
     **CANCEL**.
4. Switch to the `chem_enum_cores` table view, then back to the Markush Enumerator
   view.
   * Expected result: the **Cores** / **R-Groups** areas and the **Enumerate** button
     are still shown (the view is not empty).

## Cleanup

- Close the tables `BddMarkushCartesian` and `BddMarkushZip` without saving.

## Automation notes

- The core cards, R-group lists and the status line are the package's own DOM
  without stable `name=` anchors; the steps need package-specific readings for the
  number of cores, the number of substituents per R-number, and the status text.
- Block C step 2: the app opens pre-filled from the user's latest enumeration (here,
  the Zip run of Block B), so its content is not asserted beyond the areas being
  present.

---
{
  "order": 4,
  "datasets": ["System:AppData/SequenceTranslator/tests/chem_enum_cores.csv", "System:AppData/SequenceTranslator/tests/chem_enum_rgroups.csv"]
}
