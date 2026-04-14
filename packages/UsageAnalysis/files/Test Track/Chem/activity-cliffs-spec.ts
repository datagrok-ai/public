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

test('Chem: Activity Cliffs', async ({page}) => {
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

  await softStep('Open Activity Cliffs dialog and run', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]');
      chemMenu!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const menuItems = document.querySelectorAll('.d4-menu-item-label');
      const ac = Array.from(menuItems).find(m => m.textContent!.trim() === 'Activity Cliffs...');
      ac!.closest('.d4-menu-item')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(35000);
    const viewers = await page.evaluate(() => Array.from(grok.shell.tv.viewers).map(v => v.type));
    expect(viewers).toContain('Scatter plot');
  });

  await softStep('Click cliffs count link', async () => {
    const hasBtn = await page.evaluate(() => !!document.querySelector('[name*="button-"][name*="cliffs"]'));
    if (hasBtn) {
      await page.evaluate(async () => {
        const btn = document.querySelector('[name*="button-"][name*="cliffs"]') as HTMLElement;
        btn.click();
        await new Promise(r => setTimeout(r, 3000));
      });
      const viewerCount = await page.evaluate(() => Array.from(grok.shell.tv.viewers).length);
      expect(viewerCount).toBeGreaterThanOrEqual(3);
    }
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
