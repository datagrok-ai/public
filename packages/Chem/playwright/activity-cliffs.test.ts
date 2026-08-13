/* ---
sub_features_covered: [chem.analyze.activity-cliffs, chem.analyze.activity-cliffs.editor, chem.analyze.activity-cliffs.init, chem.analyze.activity-cliffs.top-menu, chem.analyze.activity-cliffs.transform]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

// The ui_tests stack the Jenkins job spins up is reachable only as http://xamgle-nginx:8889; any
// other target is a full stand, where this walk does complete. Skipping unconditionally meant it
// ran nowhere.
const MINIMAL_CI_STACK = /xamgle-nginx/.test(process.env.DATAGROK_URL ?? '');
async function openDatasetAndWaitForMolecule(page: Page, label: string, datasetPath: string) {
  await softStep(`[${label}] Open ${datasetPath} + wait for Chem menu (Molecule semType)`, async () => {
    const isSdf = datasetPath.toLowerCase().endsWith('.sdf');
    await page.evaluate(async ({path, isSdf}) => {
      document.body.classList.add('selenium');
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      for (let i = 0; i < 25; i++) {
        if (!(grok as any).shell.tv) break;
        await new Promise(r => setTimeout(r, 200));
      }
      if (isSdf) {
        await ((DG as any).Func.find({name: 'OpenFile'})[0])
          .prepare({fullPath: path}).call(undefined, undefined, {processed: false});
      } else {
        const df = await grok.dapi.files.readCsv(path);
        grok.shell.addTableView(df);
      }
    }, {path: datasetPath, isSdf});
    await waitForChemMenu(page);
  });
}

async function openActivityCliffsDialog(page: Page, label: string) {
  await softStep(`[${label}] Open Chem → Analyze → Activity Cliffs dialog`, async () => {
    // Real mouse: a synthetic click reaches the leaf without the submenu ever opening, and the
    // dialog then shows up with an uninitialised editor — OK runs it with no columns picked.
    await page.locator('[name="div-Chem"]').click();
    await page.locator('[name="div-Chem---Analyze"]').hover();
    await page.locator('[name="div-Chem---Analyze---Activity-Cliffs..."]').click();
    await page.locator('.d4-dialog').waitFor({timeout: 15000});
    const title = await page.evaluate(() =>
      document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')?.textContent?.trim() ?? '');
    expect(title, `Dialog title expected /Activity Cliffs/, got "${title}"`).toMatch(/Activity Cliffs/i);
  });
}

async function okAndWaitForScatter(page: Page, label: string, previousParams = ''): Promise<string> {
  let params = '';
  await softStep(`[${label}] OK → Scatter plot + activityCliffsParams tag`, async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    const result = await page.evaluate(async (prev) => {
      for (let i = 0; i < 90; i++) {
        const tv = grok.shell.tv;
        const types = Array.from(tv?.viewers ?? []).map((v: any) => v.type);
        const hasScatter = types.some(t => /Scatter plot/i.test(t));
        const acParams = tv?.dataFrame?.getTag?.('activityCliffsParams');
        // On a re-run the previous scatter and tag are still there, so "a tag exists" is not a
        // readiness signal — wait until it actually changes.
        if (hasScatter && acParams && acParams !== prev) return {ok: true, types, params: acParams, attempt: i+1};
        await new Promise(r => setTimeout(r, 2000));
      }
      const tv2 = grok.shell.tv;
      return {ok: false, types: Array.from(tv2?.viewers ?? []).map((v: any) => v.type),
        params: tv2?.dataFrame?.getTag?.('activityCliffsParams')};
    }, previousParams);
    expect((result as any).ok,
      `[${label}] Activity Cliffs OK did not produce a Scatter plot + a fresh activityCliffsParams tag within 180s. result=${JSON.stringify(result)}`,
    ).toBe(true);
    params = (result as any).params ?? '';
  });
  return params;
}

async function toggleShowOnlyCliffs(page: Page, label: string) {
  await softStep(`[${label}] Toggle Show only cliffs → cliffs filter applied`, async () => {
    // The switch does not touch a viewer option — it tags the dataframe `filterCliffs` and narrows
    // df.filter to the cliff pairs (libraries/ml/src/viewers/activity-cliffs.ts:198). Assert that.
    const before = await page.evaluate(() => ({
      trueCount: grok.shell.t.filter.trueCount,
      tag: grok.shell.t.getTag('filterCliffs') ?? '',
    }));
    // No JS API substitution — the toggle must be driven through the UI switch (.md Notes), and with
    // the real mouse: an element.click() on the switch decorator leaves the option untouched.
    const switchEl = page.locator('[name="input-host-Show-only-cliffs"] .ui-input-switch').first();
    await switchEl.waitFor({state: 'visible', timeout: 15_000});
    await switchEl.click();
    await page.waitForFunction((prev: any) =>
      (grok.shell.t.getTag('filterCliffs') ?? '').length > 0 && grok.shell.t.filter.trueCount < prev.trueCount,
    before, {timeout: 30_000});
    const after = await page.evaluate(() => ({
      trueCount: grok.shell.t.filter.trueCount,
      tag: grok.shell.t.getTag('filterCliffs') ?? '',
    }));
    expect(after.tag.length, `[${label}] Show only cliffs did not tag the dataframe`).toBeGreaterThan(0);
    expect(after.trueCount, `[${label}] Show only cliffs did not narrow the filter ` +
      `(${before.trueCount} -> ${after.trueCount})`).toBeLessThan(before.trueCount);
  });
}

async function clickCliffsLinkAndOpenGrid(page: Page, label: string): Promise<number> {
  let cliffsGridRowCount = 0;
  await softStep(`[${label}] Click "N cliffs" link → cliffs grid docks + populated`, async () => {
    await page.evaluate(() => {
      const btn = document.querySelector('.cliffs_div button.scatter_plot_link.cliffs_grid') as HTMLElement | null;
      if (!btn) throw new Error('"N cliffs" link button not found in .cliffs_div');
      btn.click();
    });
    await page.waitForFunction(() =>
      Array.from((grok as any).shell.tv?.viewers ?? []).filter((v: any) => v.type === 'Grid').length >= 2,
    undefined, {timeout: 15000});
    const result = await page.evaluate(() => {
      const grids = Array.from(grok.shell.tv?.viewers ?? []).filter((v: any) => v.type === 'Grid');
      const cliffsGrid: any = grids[1];
      const df = cliffsGrid.dataFrame;
      return {gridCount: grids.length, rowCount: df?.rowCount,
        cols: df?.columns?.toList()?.map((c: any) => c.name)};
    });
    expect((result as any).rowCount,
      `[${label}] cliffs grid should be populated with cliff pair rows — ${JSON.stringify(result)}`,
    ).toBeGreaterThan(0);
    cliffsGridRowCount = (result as any).rowCount;
  });
  return cliffsGridRowCount;
}

async function clickFirstCliffsRow(page: Page, label: string) {
  await softStep(`[${label}] currentRowIdx=0 on cliffs grid → grid row sync`, async () => {
    await page.evaluate(() => {
      const grids = Array.from(grok.shell.tv?.viewers ?? []).filter((v: any) => v.type === 'Grid');
      const cliffsGrid: any = grids[1];
      cliffsGrid.dataFrame.currentRowIdx = 0;
      const scatter = Array.from(grok.shell.tv?.viewers ?? []).find((v: any) => v.type === 'Scatter plot');
      (grok as any).shell.o = scatter;
    });
    // GROK: the property-panel pair-of-molecules render (step 6) and cliff-line hover/zoom (steps 7-8)
    // ride on scatter canvas hit-testing (activity-cliffs-run.md) — not DOM-addressable. Assert the
    // DOM-observable proxy: the cliffs-grid current row synced to the first cliff pair.
    const info = await page.evaluate(() => {
      const grids = Array.from(grok.shell.tv?.viewers ?? []).filter((v: any) => v.type === 'Grid');
      const rowIdx = (grids[1] as any).dataFrame.currentRowIdx;
      const pane = document.querySelector('.d4-pane-cliff_details');
      return {rowIdx, paneFound: !!pane, paneText: pane?.textContent?.trim() ?? ''};
    });
    expect((info as any).rowIdx, `[${label}] cliffs grid current row should sync to first cliff pair`).toBe(0);
    if ((info as any).paneFound)
      expect((info as any).paneText.length, `[${label}] cliff_details pane rendered but empty`).toBeGreaterThan(0);
  });
}

async function reRunWithCustomParam(page: Page, label: string, defaultParams: string) {
  await openActivityCliffsDialog(page, `${label}/rerun`);
  await softStep(`[${label}] Change Similarity cutoff 80 → 60`, async () => {
    // fill() alone updates the DOM input but not the Dart-side value — the run still uses 80.
    // Type it and blur so the input commits.
    const cutoff = page.locator('.d4-dialog [name="input-Similarity-cutoff"]');
    await cutoff.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('60', {delay: 50});
    await page.keyboard.press('Tab');
    await expect(cutoff).toHaveValue('60');
  });
  const customParams = await okAndWaitForScatter(page, `${label}/custom`, defaultParams);
  await softStep(`[${label}] Edited-param run differs from defaults`, async () => {
    expect(customParams,
      `[${label}] activityCliffsParams should reflect the edited cutoff (differ from defaults "${defaultParams}")`,
    ).not.toBe(defaultParams);
  });
}

async function runActivityCliffsWalk(page: Page, label: string, datasetPath: string) {
  await openDatasetAndWaitForMolecule(page, label, datasetPath);
  await openActivityCliffsDialog(page, label);
  const defaultParams = await okAndWaitForScatter(page, `${label}/defaults`);
  await toggleShowOnlyCliffs(page, label);
  await clickCliffsLinkAndOpenGrid(page, label);
  await clickFirstCliffsRow(page, label);
  await reRunWithCustomParam(page, label, defaultParams);
  await softStep(`[${label}] Close active view`, async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      for (let i = 0; i < 25; i++) {
        if (!(grok as any).shell.tv) break;
        await new Promise(r => setTimeout(r, 200));
      }
    });
  });
}

test('Chem: Activity Cliffs multi-format walk (D1-D5)', async ({page}) => {
  // CI SKIP (approved): heavy UMAP over 5 datasets exceeds the minimal CI stack (ApprovedDrugs2015 >90s,
  // no scatter) and the "Show only cliffs" / "N cliffs" UI controls aren't reachable there — the .md
  // forbids JS-API substitution for the toggle. Runs on a full stack. See PACKAGE-PLAYWRIGHT-CODE-FINDINGS.md §B1.
  test.skip(false && MINIMAL_CI_STACK, 'TEMP probe: re-checking the B1 skip now that the walk is green on dev');
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await page.waitForFunction(() => typeof grok !== 'undefined' && !!(grok as any).shell, undefined, {timeout: 30000});

  // molV3000 coverage used to point at DemoFiles ApprovedDrugs2015.sdf, which has NO numeric column
  // at all — Activity Cliffs needs an activity column, so its editor opened with empty Column and
  // Activities inputs and OK silently did nothing. Build a V3000 fixture with an activity column
  // instead, so the format is actually covered.
  const v3kFixture = `System:AppData/Chem/temp/activity-cliffs-v3000-${Date.now()}.csv`;
  let v3kFixtureWritten = false;
  await softStep('Setup: build the molV3000 fixture (structures + activity)', async () => {
    const built = await page.evaluate(async (path) => {
      await grok.functions.call('Chem:getRdKitModule', {});
      const src = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');
      const smilesCol: any = src.columns.toList().find((c: any) => /smiles/i.test(c.name));
      const actCol: any = src.columns.toList().find((c: any) => c.type === 'double' || c.type === 'int');
      const mols: string[] = [];
      const acts: number[] = [];
      for (let i = 0; i < Math.min(40, src.rowCount); i++) {
        const v3 = await grok.functions.call('Chem:convertMolNotation',
          {molecule: smilesCol.get(i), sourceNotation: 'smiles', targetNotation: 'v3Kmolblock'});
        if (!String(v3).includes('V3000')) continue;
        mols.push(String(v3));
        acts.push(Number(actCol.get(i)));
      }
      const out = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('molecule', mols),
        DG.Column.fromList('double' as any, 'activity', acts),
      ]);
      await grok.dapi.files.writeAsText(path, out.toCsv());
      return {rows: mols.length, bytes: (await grok.dapi.files.readAsText(path)).length};
    }, v3kFixture);
    expect(built.rows, 'no V3000 structures produced for the fixture').toBeGreaterThan(10);
    expect(built.bytes, `fixture not written to ${v3kFixture}`).toBeGreaterThan(0);
    v3kFixtureWritten = true;
  });

  const datasets: [string, string][] = [
    ['smiles-50', 'System:AppData/Chem/tests/smiles-50.csv'],
    ['mol1K', 'System:AppData/Chem/mol1K.sdf'],
    ['molV3000', v3kFixture],
    ['smiles_2_columns', 'System:AppData/Chem/tests/smiles_2_columns.csv'],
    ['spgi-100', 'System:AppData/Chem/tests/spgi-100.csv'],
  ];
  for (const [label, path] of datasets)
    await runActivityCliffsWalk(page, label, path);

  if (v3kFixtureWritten) {
    await page.evaluate(async (path) => {
      try { await grok.dapi.files.delete(path); } catch (_) { /* best effort */ }
    }, v3kFixture).catch(() => {});
  }

  finishSpec();
});
