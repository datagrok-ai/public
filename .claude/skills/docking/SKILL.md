---
name: docking
description: Dock small molecules to a macromolecule in Datagrok via the AutoDock GPU package
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# docking

## When to use

You have a Datagrok DataFrame of ligands (SMILES, molblock, etc.) and want
binding poses + binding-energy scores against a protein target — either
interactively from the **Chem | Docking | AutoDock...** dialog, or
programmatically from a package/script via the JS API.

## Prerequisites

- Docking package installed on the Datagrok server. Without it,
 `getAutoDockService` throws `Package 'Docking' must be installed for
 AutoDock service.` (knowledge `DG-FACT-246`).
- A **target folder** under `System:AppData/Docking/targets/<name>/`
 containing exactly two files: a `.gpf` AutoDock grid-parameter file and
 a receptor in `.pdbqt` (preferred) or `.pdb`. The folder name surfaces
 verbatim as the **Target** choice in the dialog (knowledge
 `DG-FACT-242`). Public reference targets:
 `packages/Docking/files/targets/{BACE1,kras,nmdar1}/`.
- A DataFrame column of small molecules. Cell values may be SMILES,
 molblock, or v3K molblock — the runner converts via
 `Chem:convertMolNotation` when `semType === 'Molecule'`.

## Steps

1. **Add (or pick) a target.** Drop your `<name>.pdbqt` + `<name>.gpf` into `System:AppData/Docking/targets/<name>/`. Generate them with AutoDockTools (`prepare_receptor4.py`, `prepare_gpf4.py`). The `.gpf`'s `ligand_types` line must cover every element in your ligand set — use `A C HD N NA OA SA CL` for broad coverage (see DG-FACT-243).
 ```text
 System:AppData/Docking/targets/
 BACE1/
 BACE1.pdbqt
 BACE1.gpf # ligand_types A C HD N NA OA SA CL
 ```

2. **Load ligands as a DataFrame.** Any path that yields a column with
 `semType === 'Molecule'` works.
 ```typescript
 import * as grok from 'datagrok-api/grok';
 import * as DG from 'datagrok-api/dg';

 const t = await grok.data.loadTable(
 'https://github.com/datagrok-ai/public/raw/master/'
 + 'packages/Docking/files/demo_files/demo_dataset.csv');
 await grok.data.detectSemanticTypes(t); // sets semType on SMILES col
 grok.shell.addTableView(t);
 ```
 Expected: `t.col('SMILES').semType === DG.SEMTYPE.MOLECULE`. Without
 the molecule semType the runner treats the cell as a raw string and
 skips the molblock conversion.

3. **Run docking from the UI.** Open **Chem | Docking | AutoDock...** (article shortens to "Chem > Autodock"). Fields:
 - **Ligands** — column with `semType: Molecule`.
 - **Target** — folder name from step 1.
 - **Poses** — conformations per ligand (dialog default 10; API default 30 — see DG-FACT-245).

 Expected: grid gains a `pose` column (`SEMTYPE.MOLECULE3D`) and `binding energy` (numeric, kcal/mol) — see DG-FACT-244. First run is ~1 min; same `(receptor, gpf, ligand)` tuple is cached server-side with monthly invalidation (see DG-FACT-247).

4. **Run docking from code.** The same dialog is reachable as the
 exported function `Docking:getAutodockResults` (`vectorFunc: true`,
 `action: join(table)`).
 ```typescript
 const result: DG.DataFrame = await grok.functions.call(
 'Docking:getAutodockResults', {
 table: t,
 ligands: t.col('SMILES'),
 target: 'BACE1',
 poses: 10,
 });
 // result has the input rows + 'pose' + 'binding energy' columns.
 ```
 Expected: `result.col('pose')!.semType === DG.SEMTYPE.MOLECULE3D`
 and `result.col('binding energy')` is numeric. Empty `result` (no
 columns) means the receptor/gpf wasn't readable — see step 1.

5. **Inspect a pose (Molstar zoom).** Click any cell in the `pose`
 column. The **AutoDock** widget opens in the right-hand context
 panel — embeds a Biostructure (Molstar) viewer auto-zoomed to the
 binding pocket and shows the AutoDock property table for that row.
 The widget is registered as
 `condition: Docking:isApplicableAutodock(molecule)` against
 `semType: Molecule3D` (`packages/Docking/src/package.ts:178-186`).

6. **Add extra properties / download pose.** In the AutoDock panel,
 click the `+` next to a property (intermolecular, electrostatic,
 torsional free, …) to materialize it as a column on the source
 DataFrame. To export a single pose: right-click the cell →
 **Download** → **as PDB** or **as CIF** (handlers wired by the
 BiostructureViewer package).

## Common failure modes

- **`Missing.gpf or.pdbqt/.pdb file in the target folder.`** — every target folder needs both files, extensions matching exactly (see DG-FACT-242).
- **Per-ligand failure for chlorinated compounds.** Default GPF omits `CL`. Ship a per-target `.gpf` whose `ligand_types` includes it (see DG-FACT-243).
- **Article says "Chem > Autodock", menu is empty.** Actual path is **Chem | Docking | AutoDock...**.
- **Field labelled "Conformations" in the article doesn't exist.** The dialog input is **Poses** (parameter name `poses`).
- **`Package 'Docking' must be installed for AutoDock service.`** — install/enable the Docking package (see DG-FACT-246).
- **First run hangs ~60 s with no progress.** Expected — fresh-target grid computation; don't cancel (see DG-FACT-247).
- **`getAutodockResults` returns an empty DataFrame.** Receptor parse failure or `docking-autodock` container not started (see DG-FACT-246).

## See also

- Source articles:
 - `help/develop/domains/chem/docking.md`
- Knowledge:
 - `docs/_internal/knowledge/knowledge-graph.md` — facts
 `DG-FACT-242`..`DG-FACT-247`.
- Reference packages:
 - `<public_repo>/packages/Docking/src/package.ts` — `runAutodock`,
 `getAutodockResults`, `autodockWidget`, `prepareAutoDockData`.
 - `<public_repo>/packages/Docking/src/apps/auto-dock-app.ts` —
 `AutoDockApp.init` and the per-ligand `runAutoDock` loop.
 - `<public_repo>/packages/Docking/files/targets/BACE1/` — minimal
 target folder (gpf + pdbqt) to copy when authoring a new target.
- Related skills:
 - `cheminformatics` (build the ligand DataFrame, set semType, run
 substructure / similarity filters before docking).
 - `access-data` (load the ligand CSV/SDF in the first place).
