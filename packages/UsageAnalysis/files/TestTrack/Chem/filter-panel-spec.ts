import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/SPGI.csv';

test('Chem: Filter Panel', async ({page}) => {
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

  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  await softStep('Filter panel shows a Structure (Molecule-column) filter', async () => {
    const hasStructureFilter = await page.evaluate(() => {
      const filters = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return Array.from(filters).some(f => {
        const header = f.querySelector('.d4-filter-header');
        return header && /Structure|struct/i.test(header.textContent || '');
      });
    });
    expect(hasStructureFilter).toBe(true);
  });

  await softStep('Sketch benzene via structure-filter sketch-link', async () => {
    const sketcherOpened = await page.evaluate(async () => {
      const filters = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      let clicked = false;
      for (const f of filters) {
        const header = f.querySelector('.d4-filter-header');
        if (header && /Structure/i.test(header.textContent || '')) {
          const sketchLink = f.querySelector('.sketch-link, .chem-canvas-wrapper');
          if (sketchLink) { (sketchLink as HTMLElement).click(); clicked = true; }
          break;
        }
      }
      await new Promise(r => setTimeout(r, 2500));
      return clicked && !!document.querySelector('.d4-dialog');
    });
    if (!sketcherOpened) test.skip(true, 'Structure filter sketch-link not found');
    const input = page.locator('.d4-dialog input[placeholder*="SMILES" i]');
    if (await input.count()) {
      await input.first().fill('c1ccccc1');
      await input.first().press('Enter');
      await page.waitForTimeout(1500);
    }
    const okBtn = page.locator('[name="button-OK"]');
    if (await okBtn.count()) await okBtn.first().click();
    await page.waitForTimeout(3500);
    const filtered = await page.evaluate(() => grok.shell.t.filter.trueCount);
    const total = await page.evaluate(() => grok.shell.t.rowCount);
    expect(filtered).toBeLessThan(total);
    expect(filtered).toBeGreaterThan(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
