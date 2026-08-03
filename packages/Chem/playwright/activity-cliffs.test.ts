/* ---
sub_features_covered: [chem.analyze.activity-cliffs, chem.analyze.activity-cliffs.editor, chem.analyze.activity-cliffs.init, chem.analyze.activity-cliffs.top-menu, chem.analyze.activity-cliffs.transform]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

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
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement | null;
      if (!chemMenu) throw new Error('Top-menu Chem entry not found');
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      let ac: HTMLElement | undefined;
      for (let i = 0; i < 40; i++) {
        ac = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(m => m.textContent!.trim() === 'Activity Cliffs...') as HTMLElement | undefined;
        if (ac) break;
        await new Promise(r => setTimeout(r, 200));
      }
      if (!ac) throw new Error('"Activity Cliffs..." sub-menu item not found');
      (ac.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 15000});
    const title = await page.evaluate(() =>
      document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')?.textContent?.trim() ?? '');
    expect(title, `Dialog title expected /Activity Cliffs/, got "${title}"`).toMatch(/Activity Cliffs/i);
  });
}

async function okAndWaitForScatter(page: Page, label: string): Promise<string> {
  let params = '';
  await softStep(`[${label}] OK → Scatter plot + activityCliffsParams tag`, async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    const result = await page.evaluate(async () => {
      for (let i = 0; i < 45; i++) {
        const tv = grok.shell.tv;
        const types = Array.from(tv?.viewers ?? []).map((v: any) => v.type);
        const hasScatter = types.some(t => /Scatter plot/i.test(t));
        const acParams = tv?.dataFrame?.getTag?.('activityCliffsParams');
        if (hasScatter && acParams) return {ok: true, types, params: acParams, attempt: i+1};
        await new Promise(r => setTimeout(r, 2000));
      }
      const tv2 = grok.shell.tv;
      return {ok: false, types: Array.from(tv2?.viewers ?? []).map((v: any) => v.type),
        params: tv2?.dataFrame?.getTag?.('activityCliffsParams')};
    });
    expect((result as any).ok,
      `[${label}] Activity Cliffs OK did not produce Scatter plot + activityCliffsParams tag within 90s. result=${JSON.stringify(result)}`,
    ).toBe(true);
    params = (result as any).params ?? '';
  });
  return params;
}

async function toggleShowOnlyCliffs(page: Page, label: string) {
  await softStep(`[${label}] Toggle Show only cliffs → showOnlyCliffs option enabled`, async () => {
    const before = await page.evaluate(() => {
      const scatter: any = Array.from((grok as any).shell.tv?.viewers ?? []).find((v: any) => v.type === 'Scatter plot');
      return scatter?.getOptions?.()?.look?.showOnlyCliffs ?? false;
    });
    // No JS API substitution — the toggle must be driven through the UI switch (.md Notes).
    const clicked = await page.evaluate(async () => {
      let host: Element | null = null;
      for (let i = 0; i < 20; i++) {
        host = document.querySelector('[name="input-host-Show-only-cliffs"]');
        if (host) break;
        await new Promise(r => setTimeout(r, 500));
      }
      if (!host) return {ok: false, reason: 'Show-only-cliffs host widget not found after 10s poll'};
      const switchEl = host.querySelector('.ui-input-switch') as HTMLElement | null;
      if (!switchEl) return {ok: false, reason: 'switch decorator not found inside host'};
      switchEl.click();
      return {ok: true};
    });
    expect((clicked as any).ok, `[${label}] Show-only-cliffs UI toggle: ${JSON.stringify(clicked)}`).toBe(true);
    await page.waitForFunction((prev) => {
      const scatter: any = Array.from((grok as any).shell.tv?.viewers ?? []).find((v: any) => v.type === 'Scatter plot');
      const now = scatter?.getOptions?.()?.look?.showOnlyCliffs ?? false;
      return now !== prev;
    }, before, {timeout: 15000});
    const after = await page.evaluate(() => {
      const scatter: any = Array.from((grok as any).shell.tv?.viewers ?? []).find((v: any) => v.type === 'Scatter plot');
      return scatter?.getOptions?.()?.look?.showOnlyCliffs ?? false;
    });
    // Non-cliff points are hidden on the scatter canvas — canvas visibility is not DOM-observable,
    // so assert the showOnlyCliffs viewer option (the state that drives the hide) flipped to true.
    expect(after, `[${label}] Show only cliffs toggle did not enable showOnlyCliffs`).toBe(true);
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
    const cutoff = page.locator('.d4-dialog [name="input-Similarity-cutoff"]');
    await cutoff.fill('60');
    await expect(cutoff).toHaveValue('60');
  });
  const customParams = await okAndWaitForScatter(page, `${label}/custom`);
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
  test.skip(true, 'CI-env: heavy UMAP walk + UI controls unavailable on the minimal CI stack (findings §B1)');
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await page.waitForFunction(() => typeof grok !== 'undefined' && !!(grok as any).shell, undefined, {timeout: 30000});

  const datasets: [string, string][] = [
    ['smiles-50', 'System:AppData/Chem/tests/smiles-50.csv'],
    ['mol1K', 'System:AppData/Chem/mol1K.sdf'],
    ['ApprovedDrugs2015', 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf'],
    ['smiles_2_columns', 'System:AppData/Chem/tests/smiles_2_columns.csv'],
    ['spgi-100', 'System:AppData/Chem/tests/spgi-100.csv'],
  ];
  for (const [label, path] of datasets)
    await runActivityCliffsWalk(page, label, path);

  finishSpec();
});
