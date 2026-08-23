/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openChemMenuItem} from '../helpers/chem';

declare const grok: any;
declare const DG: any;

// Selector recon-notes (class-2: derived from the registered function signatures in
// public/packages/Chem/src/package.g.ts, read 2026-08-22; added to chem.md in the same change.
// Source-derived, not live-observed — every name below follows the documented
// caption→`[name="input-<Caption>"]` / title→`[name="dialog-<Title>"]` convention):
//   [name="dialog-BitBIRCH-Clustering"] — the top-menu dialog of `Chem | Calculate | BitBIRCH
//     Clustering...` (package.g.ts:345-354). Inputs: [name="input-host-Molecules"] (column, no
//     caption — host named after the `molecules` input), [name="input-host-Threshold"] (double,
//     registered default 0.55), [name="input-host-Fingerprint-type"] (choice, default 'Morgan').
//   [name="dialog-Cluster-MCS"] — top-menu dialog of `Chem | Calculate | Cluster MCS...`
//     (package.g.ts:356-364). Inputs: [name="input-host-Molecules"] (caption Molecules) and
//     [name="input-host-Cluster"] (caption Cluster, `nullable: false` — required, so the platform
//     pre-selects the first column passing the input's filter). That filter is `type: 'categorical'`
//     (package.ts:833), which resolves to `Column.isCategorical` (column_filter.dart:56) — true only
//     for string and bool columns (string_column.dart:23, bool_column.dart:13; int_column.dart:17
//     defaults it to false and nothing sets it). The BitBIRCH cluster-id column of S1 is an int
//     column (bitbirch-clustering.ts:37,58), so it is NOT offered here. In smiles.csv the only
//     categorical column is `canonical_smiles` itself — see the S3.4 note on what that implies.
//   [name="dialog-Similarity-Matrix"] — top-menu dialog of `Chem | Calculate | Similarity
//     Matrix...` (package.g.ts:235-243). Inputs: [name="input-host-Molecules"],
//     [name="input-host-Symbols"], [name="input-host-Fingerprint-type"].
//   Column value inside a column-selector host is read from `.d4-column-selector-column`, the same
//     shape chem.md documents for the Descriptors / Chemical Properties dialogs.
//   Result-table name of Similarity Matrix is `<molecule column name> similarity matrix` — it is
//     built from the COLUMN name, not the table name (package.ts:603).

test.use(specTestOptions);

const MOL_COL = 'canonical_smiles';

async function columnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.t.columns.names());
}

async function waitForNewColumn(page: Page, before: string[], match: string, capMs: number): Promise<string | null> {
  await page.waitForFunction(({before, match}) => {
    const re = new RegExp(match, 'i');
    return grok.shell.t.columns.names().some((n: string) => !before.includes(n) && re.test(n));
  }, {before, match}, {timeout: capMs}).catch(() => {});
  const after = await columnNames(page);
  const re = new RegExp(match, 'i');
  return after.find((n) => !before.includes(n) && re.test(n)) ?? null;
}

test('Chem: Calculate — BitBIRCH Clustering, Cluster MCS, Similarity Matrix', async ({page}) => {
  test.setTimeout(600_000);

  const consoleErrors: string[] = [];
  const AMBIENT = [/favicon/i, /ResizeObserver loop/i, /Permissions policy violation/i,
    /Unable to find element in cloned iframe/i];
  const isAmbient = (text: string) => AMBIENT.some((re) => re.test(text));
  page.on('console', (msg) => { if (msg.type() === 'error' && !isAmbient(msg.text())) consoleErrors.push(msg.text()); });
  page.on('pageerror', (e) => { if (!isAmbient(String(e))) consoleErrors.push(String(e)); });
  const expectNoConsoleErrors = (where: string) => {
    const seen = consoleErrors.splice(0, consoleErrors.length);
    expect(seen, `JavaScript console errors fired during ${where}`).toEqual([]);
  };

  await loginToDatagrok(page);

  let baseline = 0;
  let beforeBitbirch: string[] = [];
  let beforeMcs: string[] = [];
  let clusterColName: string | null = null;

  await softStep('Setup: open smiles.csv, wait for Molecule semType + Chem menu; record baseline row count', async () => {
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
      const detected = new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 5000);
      });
      grok.shell.addTableView(df);
      await detected;
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await page.locator('[name="viewer-Grid"] canvas').first().waitFor({timeout: 30000, state: 'attached'});
    await waitForChemMenu(page);
    await waitForMolecule(page);
    baseline = await page.evaluate(() => grok.shell.t.rowCount);
    expect(baseline, 'smiles.csv baseline row count must be established and positive').toBeGreaterThan(0);
    const molTyped = await page.evaluate((c) => grok.shell.t.col(c)?.semType === 'Molecule', MOL_COL);
    expect(molTyped, `the ${MOL_COL} column must be typed as Molecule before any clustering runs`).toBe(true);
  });

  await softStep('S1.1-3: Chem → Calculate → BitBIRCH Clustering → dialog opens, accept defaults, OK', async () => {
    await openChemMenuItem(page, 'BitBIRCH Clustering...', {delayMs: 800});
    await page.locator('[name="dialog-BitBIRCH-Clustering"]').waitFor({timeout: 15000});
    const molHost = await page.locator('[name="dialog-BitBIRCH-Clustering"] [name="input-host-Molecules"] .d4-column-selector-column').textContent();
    expect(molHost?.trim(), 'the molecule column selector must auto-detect the Molecule column').toBe(MOL_COL);
    await expect(page.locator('[name="dialog-BitBIRCH-Clustering"] [name="input-host-Threshold"]'),
      'the BitBIRCH dialog must expose a Threshold control — the clustering cut-off is the one parameter the scenario ' +
      'accepts by default, and a dialog without it is not the BitBIRCH editor').toBeVisible();
    beforeBitbirch = await columnNames(page);
    await page.locator('[name="dialog-BitBIRCH-Clustering"] [name="button-OK"]').click();
  });

  await softStep('S1.4: a cluster-id column is appended; distinct-value count is > 1 and < row count; row count unchanged; no console errors', async () => {
    const col = await waitForNewColumn(page, beforeBitbirch, 'cluster|bitbirch', 90000);
    expect(col, 'a BitBIRCH cluster-id column must appear within 90s of the OK click').toBeTruthy();
    const probe = await page.evaluate(({baseline, name}) => {
      const t = grok.shell.t;
      const c = name ? t.col(name) : null;
      const distinct = c ? new Set(Array.from({length: t.rowCount}, (_: any, i: number) => c.get(i))).size : 0;
      return {name, distinct, rowCount: t.rowCount, baseline};
    }, {baseline, name: col});
    expect(probe.name, 'a BitBIRCH cluster-id column must be appended').toBeTruthy();
    expect(probe.distinct,
      'the cluster-id column must hold more than one distinct value — a single cluster means BitBIRCH did not partition').toBeGreaterThan(1);
    expect(probe.distinct,
      'the distinct-cluster count must be strictly less than the row count — one cluster per molecule is a trivial non-clustering').toBeLessThan(probe.rowCount);
    expect(probe.rowCount, 'the grid row count must be unchanged after the cluster-id column-append').toBe(probe.baseline);
    expectNoConsoleErrors('the BitBIRCH dialog open, OK click, and column-append');
  });

  await softStep('S3.1-3: Chem → Calculate → Cluster MCS → dialog opens, Molecule detected, OK', async () => {
    await openChemMenuItem(page, 'Cluster MCS...', {delayMs: 800});
    await page.locator('[name="dialog-Cluster-MCS"]').waitFor({timeout: 15000});
    const molHost = await page.locator('[name="dialog-Cluster-MCS"] [name="input-host-Molecules"] .d4-column-selector-column').textContent();
    expect(molHost?.trim(), 'the Cluster MCS molecule selector must show the Molecule column').toBe(MOL_COL);
    const clusterHost = await page.locator('[name="dialog-Cluster-MCS"] [name="input-host-Cluster"] .d4-column-selector-column').textContent();
    // `.d4-column-selector-column` renders `column.friendlyName` (column_combo_box.dart:126), not the
    // column identity the call is made with. The effective value is the ColumnInput's `stringValue`,
    // which is `value?.name` (column_input.dart:53) — reachable through the live dialog model.
    const clusterEffective = await page.evaluate(() => {
      const dlg = DG.Dialog.getOpenDialogs().find((d: any) => d.title === 'Cluster MCS');
      const input = dlg?.inputs.find((i: any) => i.caption === 'Cluster');
      return {value: input?.stringValue ?? null, inputType: input?.inputType ?? null};
    });
    console.log(`[S3.1] Cluster MCS dialog: Molecules="${molHost?.trim()}", Cluster rendered="${clusterHost?.trim()}", ` +
      `Cluster effective value="${clusterEffective.value}" (input type ${clusterEffective.inputType})`);
    expect(clusterHost?.trim(),
      'the Cluster MCS dialog must expose a Cluster column selector carrying a pre-selected column — the input is ' +
      'registered nullable: false, so an empty selector would block OK and the default-accepting run below could not proceed')
      .toBeTruthy();
    clusterColName = clusterEffective.value;
    beforeMcs = await columnNames(page);
    await page.locator('[name="dialog-Cluster-MCS"] [name="button-OK"]').click();
  });

  await softStep('S3.4: a Cluster MCS column is appended, semType Molecule, every row carrying a non-empty structure string; row count unchanged', async () => {
    const mcsCol = await waitForNewColumn(page, beforeMcs, 'mcs|scaffold', 120000);
    expect(mcsCol, 'a Cluster MCS column must appear within 120s of the OK click').toBeTruthy();
    const probe = await page.evaluate(({baseline, name, molCol, clusterName}) => {
      const t = grok.shell.t;
      const c = name ? t.col(name) : null;
      const src = t.col(molCol);
      const cluster = clusterName ? t.col(clusterName) : null;
      const vals = c ? Array.from({length: t.rowCount}, (_: any, i: number) => c.get(i)) : [];
      const nonEmpty = vals.filter((v: any) => typeof v === 'string' && v.trim().length > 0).length;
      const canonicalized = vals.filter((v: any, i: number) =>
        typeof v === 'string' && v.trim().length > 0 && v !== src.get(i)).length;
      const sizes = new Map<string, number>();
      if (cluster)
        for (let i = 0; i < t.rowCount; i++) {
          const k = String(cluster.get(i));
          sizes.set(k, (sizes.get(k) ?? 0) + 1);
        }
      const counts = Array.from(sizes.values());
      return {name, semType: c?.semType ?? null, nonEmpty, canonicalized, rowCount: t.rowCount, baseline,
        clusterName, clusterCount: counts.length, largestCluster: counts.length ? Math.max(...counts) : 0};
    }, {baseline, name: mcsCol, molCol: MOL_COL, clusterName: clusterColName});
    expect(probe.name, 'a Cluster MCS column must be appended').toBeTruthy();
    console.log(`[S3.4] MCS column ${probe.name}: semType=${probe.semType}, ${probe.nonEmpty}/${probe.rowCount} ` +
      `non-empty cells, ${probe.canonicalized} of them differ textually from the row's own ${MOL_COL} value; ` +
      `partition column "${probe.clusterName}" holds ${probe.clusterCount} clusters over ${probe.rowCount} rows, ` +
      `largest cluster ${probe.largestCluster} member(s)`);
    // SCOPE REDUCTION SR-01 (see the scenario's scope_reductions). Nothing asserted here is about a
    // computed maximum common substructure. The Cluster input only accepts categorical columns
    // (package.ts:833 -> Column.isCategorical), and smiles.csv's only categorical column is
    // `canonical_smiles` itself, whose 1000 values are all distinct — so every cluster is a singleton
    // and the multi-member MCS branch (rdkit-service.ts:633-648) never runs. `canonicalized` is not a
    // discriminator either: a singleton takes `mols[0]` (rdkit-service.ts:632) but every entry is then
    // re-parsed and re-emitted through `get_smiles()` (rdkit-service.ts:659-669), so a passed-through
    // molecule comes back RDKit-canonicalized and differs from the source whenever the source spelling
    // was not already canonical. It is logged as evidence of the round-trip, never asserted on.
    // What the count below does discriminate: the failure paths that leave entries empty — the
    // pre-filled '' array (rdkit-service.ts:618), the 45s worker-timeout restart (:634-645), the
    // per-worker catch (:643-645), and a get_smiles() parse failure (:667).
    // The premise of SR-01, asserted rather than merely logged: a reduction whose premise can quietly
    // stop holding is a silent false-green in the documentation layer.
    expect(probe.largestCluster,
      'the partition pre-selected by the Cluster MCS dialog must put every molecule in a cluster of one — ' +
      'this is not a claim about the product but the premise of scope reduction SR-01, which excuses this ' +
      'step from checking the value written for a multi-member cluster precisely because no multi-member ' +
      'cluster is reachable here. A larger cluster means the premise has gone stale (a fixture change, or a ' +
      'change to which columns the categorical filter offers) and SR-01 must be revisited and this step ' +
      'strengthened to check the computed common substructure — not that a defect should be filed. A value ' +
      'of 0 means the partition column was never resolved, so the premise was not measured at all').toBe(1);
    expect(probe.nonEmpty,
      'every row of the MCS column must carry a non-empty structure string — the product pre-fills the ' +
      'per-cluster result with empty strings and leaves them empty on worker timeout, worker error, or a ' +
      'failed RDKit re-parse, so any empty cell is one of those failure paths').toBe(probe.rowCount);
    // semType is set only on the success path (package.ts:860); the failure handler returns the column
    // untouched, so a null semType is the observable of a swallowed 'Cluster MCS Error'
    expect(probe.semType,
      'the appended MCS column must carry semType Molecule — the product tags it only after the scaffold ' +
      'computation returns, so an untagged column is a failed run that was swallowed into an empty column').toBe('Molecule');
    expect(probe.rowCount, 'the grid row count must be unchanged after the MCS column-append').toBe(probe.baseline);
    expectNoConsoleErrors('the Cluster MCS dialog open, OK click, and column-append');
  });

  await softStep('S4 setup: extract a 50-row Molecule subset so the pairwise matrix is bounded, activate it', async () => {
    await page.evaluate((molCol) => {
      const full = grok.shell.t;
      const smiles = full.col(molCol);
      const n = Math.min(50, full.rowCount);
      const small = DG.DataFrame.create(n);
      const scol = DG.Column.fromStrings(molCol, Array.from({length: n}, (_: any, i: number) => smiles.get(i)));
      scol.semType = 'Molecule';
      small.columns.add(scol);
      small.name = 'clustering_matrix_subset';
      grok.shell.addTableView(small);
    }, MOL_COL);
    await page.waitForFunction(() => grok.shell.t.name === 'clustering_matrix_subset' &&
      grok.shell.t.col('canonical_smiles')?.semType === 'Molecule', null, {timeout: 20000});
    await waitForChemMenu(page);
    const subsetRows = await page.evaluate(() => grok.shell.t.rowCount);
    expect(subsetRows, 'the extracted subset must have rows to build a matrix from').toBeGreaterThan(1);
  });

  await softStep('S4.1-3: Chem → Calculate → Similarity Matrix → accept defaults, OK', async () => {
    await openChemMenuItem(page, 'Similarity Matrix...', {delayMs: 800});
    await page.locator('[name="dialog-Similarity-Matrix"]').waitFor({timeout: 15000});
    const molHost = await page.locator('[name="dialog-Similarity-Matrix"] [name="input-host-Molecules"] .d4-column-selector-column').textContent();
    expect(molHost?.trim(), 'the Similarity Matrix molecule selector must show the Molecule column').toBe(MOL_COL);
    await page.locator('[name="dialog-Similarity-Matrix"] [name="button-OK"]').click();
  });

  await softStep('S4.4: a similarity-matrix table named after the molecule column is produced; the diagonal reads 1.0 (row/column alignment); at least one off-diagonal cell differs from 1.0', async () => {
    await page.waitForFunction(() => grok.shell.tables.some((t: any) => /similarity matrix/i.test(t.name)),
      null, {timeout: 120000}).catch(() => {});
    const probe = await page.evaluate(() => {
      const matrix = grok.shell.tables.find((t: any) => /similarity matrix/i.test(t.name));
      if (!matrix) return {found: false};
      const floatCols = matrix.columns.names().filter((n: string) => n !== 'symbol');
      const n = Math.min(matrix.rowCount, floatCols.length);
      let diagAllOne = true;
      let offDiagDiffers = false;
      for (let i = 0; i < n; i++) {
        const c = matrix.col(floatCols[i]);
        if (Math.abs(c.get(i) - 1.0) > 1e-6) diagAllOne = false;
        const j = (i + 1) % matrix.rowCount;
        if (Math.abs(c.get(j) - 1.0) > 1e-6) offDiagDiffers = true;
      }
      return {found: true, name: matrix.name, rowCount: matrix.rowCount, floatColCount: floatCols.length, diagAllOne, offDiagDiffers};
    });
    expect(probe.found, 'a similarity-matrix result table must be produced').toBe(true);
    console.log(`[S4.4] matrix table "${probe.name}": ${probe.rowCount} rows x ${probe.floatColCount} similarity columns`);
    expect(probe.name,
      `the matrix table must be named after the molecule column it was built from (\`${MOL_COL} similarity matrix\`) — ` +
      'any other name means the probe latched onto a table this run did not produce').toMatch(new RegExp(`^${MOL_COL} similarity matrix`));
    expect(probe.rowCount,
      'the matrix must cover more than one molecule, otherwise a 1x1 diagonal proves nothing').toBeGreaterThan(1);
    expect(probe.floatColCount,
      'the matrix must carry exactly one similarity column per molecule — otherwise column i is not molecule i and the cells read below are not the diagonal').toBe(probe.rowCount);
    // the product writes the diagonal as the literal 1 (package.ts:609), never as a computed
    // self-similarity, so this guards the row/column alignment of the result, not the metric
    expect(probe.diagAllOne,
      'cell [i][i] of the matrix must read 1.0 — the diagonal is written as a literal, so a value off 1.0 ' +
      'means column i is not molecule i and every cell read by row/column index is misaligned').toBe(true);
    expect(probe.offDiagDiffers,
      'at least one off-diagonal cell must differ from 1.0 — an all-1.0 matrix means no real pairwise comparison was computed').toBe(true);
    expectNoConsoleErrors('the Similarity Matrix run and its viewer rendering');
  });

  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
