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
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await waitForChemMenu(page);
  });
}

async function openChemicalSpaceDialog(page: Page, label: string) {
  await softStep(`[${label}] Open Chem → Analyze → Chemical Space dialog`, async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement | null;
      if (!chemMenu) throw new Error('Top-menu Chem entry not found');
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 800));
      const cs = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Chemical Space...') as HTMLElement | undefined;
      if (!cs) throw new Error('"Chemical Space..." sub-menu item not found');
      (cs.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 15000});
    const title = await page.evaluate(() =>
      document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')?.textContent?.trim() ?? '');
    expect(title, `Dialog title expected /Chemical/, got "${title}"`).toMatch(/Chemical/i);
  });
}

async function clickOkAndWaitForEmbedding(page: Page, label: string, minSuffix: number) {
  let foundSuffix = -1;
  await softStep(`[${label}] OK → embedding scatter plot + Embed_X_N (N > ${minSuffix})`, async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    const result = await page.evaluate(async (lastSuffix: number) => {
      for (let i = 0; i < 45; i++) {
        const df = grok.shell.tv?.dataFrame;
        if (!df) { await new Promise(r => setTimeout(r, 2000)); continue; }
        const embedCols = df.columns.toList().filter((c: any) => /^Embed_X_(\d+)$/.test(c.name));
        const hasScatter = Array.from(grok.shell.tv.viewers).some((v: any) => v.type === 'Scatter plot');
        if (embedCols.length > 0 && hasScatter) {
          const maxSuffix = Math.max(...embedCols.map((c: any) => parseInt(c.name.match(/^Embed_X_(\d+)$/)![1])));
          if (maxSuffix > lastSuffix) {
            // Tag is set asynchronously ~5-8s after column creation; poll up to 20s.
            let hasTag = false;
            for (let j = 0; j < 10; j++) {
              const newCol = grok.shell.tv?.dataFrame?.col(`Embed_X_${maxSuffix}`);
              if (newCol && Array.from(newCol.tags.keys()).some((k: any) => /chem-space-embedding-col/.test(k))) {
                hasTag = true;
                break;
              }
              await new Promise(r => setTimeout(r, 2000));
            }
            return {ok: true, maxSuffix, hasTag, viewerTypes: Array.from(grok.shell.tv.viewers).map((v: any) => v.type)};
          }
        }
        await new Promise(r => setTimeout(r, 2000));
      }
      const df = grok.shell.tv?.dataFrame;
      return {ok: false, dfCols: df?.columns?.toList()?.map((c: any) => c.name)?.slice(0, 30),
        viewerTypes: Array.from(grok.shell.tv?.viewers ?? []).map((v: any) => v.type)};
    }, minSuffix);
    expect((result as any).ok,
      `[${label}] Chemical Space did not produce a fresh Embed_X_N (N > ${minSuffix}) within 90s. result=${JSON.stringify(result)}`,
    ).toBe(true);
    expect((result as any).hasTag,
      `[${label}] new Embed_X column missing '.%chem-space-embedding-col' tag after 20s wait`,
    ).toBe(true);
    foundSuffix = (result as any).maxSuffix;
  });
  return foundSuffix;
}

async function runChemicalSpaceWalk(page: Page, label: string, datasetPath: string,
  customParam: 'method' | 'cluster-mcs') {

  await openDatasetAndWaitForMolecule(page, label, datasetPath);
  await openChemicalSpaceDialog(page, label);
  const defaultSuffix = await clickOkAndWaitForEmbedding(page, `${label}/defaults`, 0);
  await openChemicalSpaceDialog(page, `${label}/reopen`);

  await softStep(`[${label}] Change param: ${customParam}`, async () => {
    if (customParam === 'method')
      await page.locator('.d4-dialog [name="input-Method"]').selectOption('t-SNE');
    else
      await page.locator('.d4-dialog [name="input-Cluster-MCS"]').check();
    await page.waitForTimeout(400);
  });

  await clickOkAndWaitForEmbedding(page, `${label}/custom`, defaultSuffix);

  await softStep(`[${label}] Close active view`, async () => {
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1500);
  });
}

test('Chem: Chemical Space multi-format walk (smiles-50 / molV2000 / molV3000)', async ({page}) => {
  test.setTimeout(900_000); // 15 min for 3 × 7-step walks on cold session

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  await runChemicalSpaceWalk(page, 'D1 smiles-50', 'System:AppData/Chem/tests/smiles-50.csv', 'method');
  await runChemicalSpaceWalk(page, 'D2 molV2000', 'System:AppData/Chem/mol1K.sdf', 'method');
  await runChemicalSpaceWalk(page, 'D3 molV3000', 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf', 'cluster-mcs');

  finishSpec();
});
