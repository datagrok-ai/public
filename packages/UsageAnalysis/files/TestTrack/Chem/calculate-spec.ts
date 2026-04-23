import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/chem/smiles.csv';

test('Chem: Calculate Descriptors', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv(path);
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
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open Chem → Calculate → Descriptors dialog', async () => {
    await page.locator('[name="div-Chem"]').click();
    await page.waitForTimeout(500);
    await page.locator('.d4-menu-item-label', {hasText: /^Calculate$/}).first().hover();
    await page.waitForTimeout(500);
    await page.locator('.d4-menu-item-label', {hasText: /^Descriptors\.\.\.$/}).first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
  });

  await softStep('Accept defaults → OK → new descriptor columns appended', async () => {
    const colsBefore = await page.evaluate(() => grok.shell.t.columns.length);
    // Expand a descriptor group first so at least one is selected, then OK.
    await page.evaluate(async () => {
      // Check any tree node (e.g., first descriptor category like "RDKit descriptors") to enable OK
      const tree = document.querySelector('.d4-dialog .d4-tree-view-root');
      if (tree) {
        const firstGroup = tree.querySelector('.d4-tree-view-checkbox');
        if (firstGroup) (firstGroup as HTMLElement).click();
      }
      await new Promise(r => setTimeout(r, 500));
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(20000);
    const colsAfter = await page.evaluate(() => grok.shell.t.columns.length);
    expect(colsAfter).toBeGreaterThan(colsBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
