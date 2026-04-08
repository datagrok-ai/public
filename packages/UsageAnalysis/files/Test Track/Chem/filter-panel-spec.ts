import {test, expect} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/SPGI.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Chem: Filter Panel', async ({page}) => {
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell && document.querySelector('.d4-root'), {timeout: 30000});

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch(e) {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Open filter panel
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  await softStep('Sketch benzene in structure filter', async () => {
    // Click sketch-link on Structure filter
    await page.evaluate(async () => {
      const filterPanel = document.querySelector('[name="viewer-Filters"]');
      const filters = filterPanel!.querySelectorAll('.d4-filter');
      for (const f of filters) {
        const header = f.querySelector('.d4-filter-header');
        if (header?.textContent?.trim() === 'Structure') {
          const sketchLink = f.querySelector('.sketch-link');
          if (sketchLink) (sketchLink as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 2000));
    });

    // Type SMILES
    const input = page.locator('input[placeholder*="SMILES"]');
    await input.fill('c1ccccc1');
    await input.press('Enter');
    await page.waitForTimeout(1000);
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(3000);

    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeLessThan(3624);
    expect(filtered).toBeGreaterThan(0);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
