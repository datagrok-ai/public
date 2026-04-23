import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function openDataset(page: any, path: string) {
  await page.evaluate(async (p: string) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv(p);
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
  }, path);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

async function runRGroupsDialog(page: any) {
  await page.evaluate(async () => {
    const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
    chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    await new Promise(r => setTimeout(r, 500));
    const item = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find(m => m.textContent!.trim() === 'R-Groups Analysis...') as HTMLElement;
    (item.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('.d4-dialog').waitFor({timeout: 10000});
}

test('Chem: R-Groups Analysis', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await softStep('Open sar_small.csv and run R-Groups with MCS', async () => {
    await openDataset(page, 'System:DemoFiles/chem/sar_small.csv');
    await runRGroupsDialog(page);
    await page.evaluate(async () => {
      const mcs = Array.from(document.querySelectorAll('.d4-dialog button'))
        .find(b => b.textContent!.trim() === 'MCS') as HTMLElement;
      mcs.click();
      await new Promise(r => setTimeout(r, 8000));
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(12000);
    const result = await page.evaluate(() => ({
      viewers: Array.from(grok.shell.tv.viewers).map((v: any) => v.type),
      cols: grok.shell.t.columns.names(),
    }));
    expect(result.viewers).toContain('Trellis plot');
    expect(result.cols.some((n: string) => /^R[1-4]$/.test(n))).toBe(true);
  });

  await softStep('smiles.csv → R-Groups produces no Trellis (no R Groups found)', async () => {
    await openDataset(page, 'System:DemoFiles/chem/smiles.csv');
    await runRGroupsDialog(page);
    await page.evaluate(async () => {
      const mcs = Array.from(document.querySelectorAll('.d4-dialog button'))
        .find(b => b.textContent!.trim() === 'MCS') as HTMLElement;
      mcs.click();
      await new Promise(r => setTimeout(r, 8000));
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(10000);
    const viewers = await page.evaluate(() => Array.from(grok.shell.tv.viewers).map((v: any) => v.type));
    expect(viewers).not.toContain('Trellis plot');
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
