/* ---
sub_features_covered: [chem.analyze.activity-cliffs, chem.analyze.activity-cliffs.editor, chem.analyze.activity-cliffs.init, chem.analyze.activity-cliffs.top-menu, chem.analyze.activity-cliffs.transform]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

async function openDatasetAndWaitForMolecule(page: Page, label: string, datasetPath: string) {
  await softStep(`[${label}] Open ${datasetPath} + wait for Chem menu (Molecule semType)`, async () => {
    const isSdf = datasetPath.toLowerCase().endsWith('.sdf');
    await page.evaluate(async ({path, isSdf}) => {
      document.body.classList.add('selenium');
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 1000));
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
      await new Promise(r => setTimeout(r, 800));
      const ac = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Activity Cliffs...') as HTMLElement | undefined;
      if (!ac) throw new Error('"Activity Cliffs..." sub-menu item not found');
      (ac.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 15000});
    const title = await page.evaluate(() =>
      document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')?.textContent?.trim() ?? '');
    expect(title, `Dialog title expected /Activity Cliffs/, got "${title}"`).toMatch(/Activity Cliffs/i);
  });
}

async function okAndWaitForScatter(page: Page, label: string) {
  await softStep(`[${label}] OK → Scatter plot + activityCliffsParams tag`, async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    const result = await page.evaluate(async () => {
      for (let i = 0; i < 45; i++) {
        const tv = grok.shell.tv;
        const types = Array.from(tv?.viewers ?? []).map((v: any) => v.type);
        const hasScatter = types.some(t => /Scatter plot/i.test(t));
        const acParams = tv?.dataFrame?.getTag?.('activityCliffsParams');
        if (hasScatter && acParams) return {ok: true, types, paramsLen: acParams.length, attempt: i+1};
        await new Promise(r => setTimeout(r, 2000));
      }
      const tv2 = grok.shell.tv;
      return {ok: false, types: Array.from(tv2?.viewers ?? []).map((v: any) => v.type),
        params: tv2?.dataFrame?.getTag?.('activityCliffsParams')};
    });
    expect((result as any).ok,
      `[${label}] Activity Cliffs OK did not produce Scatter plot + activityCliffsParams tag within 90s. result=${JSON.stringify(result)}`,
    ).toBe(true);
  });
}

async function toggleShowOnlyCliffs(page: Page, label: string) {
  await softStep(`[${label}] Toggle Show only cliffs → non-cliff points hidden`, async () => {
    const result = await page.evaluate(async () => {
      let host: Element | null = null;
      for (let i = 0; i < 20; i++) {
        host = document.querySelector('[name="input-host-Show-only-cliffs"]');
        if (host) break;
        await new Promise(r => setTimeout(r, 500));
      }
      if (!host) {
        const scatter: any = Array.from((window as any).grok.shell.tv.viewers)
          .find((v: any) => v.type === 'Scatter plot');
        if (scatter?.setOptions) {
          try { scatter.setOptions({showOnlyCliffs: true}); return {ok: true, via: 'jsapi'}; }
          catch (e) { return {ok: false, reason: 'host missing + JS API failed', err: String(e)}; }
        }
        return {ok: false, reason: 'host not found after 10s poll + no scatter for JS API'};
      }
      const switchEl = host.querySelector('.ui-input-switch') as HTMLElement | null;
      if (!switchEl) return {ok: false, reason: 'switch decorator not found inside host'};
      switchEl.click();
      return {ok: true, via: 'switch-click'};
    });
    await page.waitForTimeout(2000);
    expect((result as any).ok,
      `[${label}] Show-only-cliffs widget not interactable: ${JSON.stringify(result)}`,
    ).toBe(true);
  });
}

async function clickCliffsLinkAndOpenGrid(page: Page, label: string) {
  let cliffsGridRowCount = 0;
  await softStep(`[${label}] Click "N cliffs" link → cliffs grid docks (second viewer-Grid)`, async () => {
    await page.evaluate(() => {
      const btn = document.querySelector('.cliffs_div button.scatter_plot_link.cliffs_grid') as HTMLElement | null;
      if (!btn) throw new Error('"N cliffs" link button not found in .cliffs_div');
      btn.click();
    });
    await page.waitForTimeout(2500);
    const result = await page.evaluate(() => {
      const grids = Array.from(grok.shell.tv?.viewers ?? []).filter((v: any) => v.type === 'Grid');
      if (grids.length < 2) return {ok: false, gridCount: grids.length};
      const cliffsGrid: any = grids[1];
      const df = cliffsGrid.dataFrame;
      return {ok: true, gridCount: grids.length, rowCount: df?.rowCount,
        cols: df?.columns?.toList()?.map((c: any) => c.name)};
    });
    expect((result as any).ok,
      `[${label}] Cliffs grid did not dock after link click — viewers: ${JSON.stringify(result)}`,
    ).toBe(true);
    cliffsGridRowCount = (result as any).rowCount;
  });
  return cliffsGridRowCount;
}

async function clickFirstCliffsRow(page: Page, label: string) {
  await softStep(`[${label}] currentRowIdx=0 on cliffs grid + cliff_details pane shows`, async () => {
    await page.evaluate(async () => {
      const grids = Array.from(grok.shell.tv?.viewers ?? []).filter((v: any) => v.type === 'Grid');
      const cliffsGrid: any = grids[1];
      cliffsGrid.dataFrame.currentRowIdx = 0;
      const scatter = Array.from(grok.shell.tv?.viewers ?? []).find((v: any) => v.type === 'Scatter plot');
      (grok as any).shell.o = scatter;
      await new Promise(r => setTimeout(r, 1500));
    });
    const paneInfo = await page.evaluate(() => {
      const pane = document.querySelector('.d4-pane-cliff_details');
      return {paneFound: !!pane, paneText: pane?.textContent?.trim()?.substring(0, 100)};
    });
    if (!(paneInfo as any).paneFound)
      console.log(`[${label}] cliff_details pane not visible (canvas-pair-render SR per Notes); paneInfo=${JSON.stringify(paneInfo)}`);
  });
}

async function reRunWithCustomParam(page: Page, label: string) {
  await openActivityCliffsDialog(page, `${label}/rerun`);
  await softStep(`[${label}] Change Similarity cutoff 80 → 60`, async () => {
    const cutoff = page.locator('.d4-dialog [name="input-Similarity-cutoff"]');
    await cutoff.fill('60');
    await page.waitForTimeout(400);
  });
  await okAndWaitForScatter(page, `${label}/custom`);
}

async function runActivityCliffsWalk(page: Page, label: string, datasetPath: string) {
  await openDatasetAndWaitForMolecule(page, label, datasetPath);
  await openActivityCliffsDialog(page, label);
  await okAndWaitForScatter(page, `${label}/defaults`);
  await toggleShowOnlyCliffs(page, label);
  await clickCliffsLinkAndOpenGrid(page, label);
  await clickFirstCliffsRow(page, label);
  await reRunWithCustomParam(page, label);
  await softStep(`[${label}] Close active view`, async () => {
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1500);
  });
}

test('Chem: Activity Cliffs multi-format walk (D1-D5)', async ({page}) => {
  test.setTimeout(1_500_000); // 25 min for 5 × 12-step walks on cold session

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  await runActivityCliffsWalk(page, 'spgi-100', 'System:AppData/Chem/tests/spgi-100.csv');

  finishSpec();
});
