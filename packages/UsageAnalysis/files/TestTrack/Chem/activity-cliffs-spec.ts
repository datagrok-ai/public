import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Chem: Activity Cliffs', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
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
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open Chem → Analyze → Activity Cliffs → dialog opens', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const ac = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Activity Cliffs...') as HTMLElement;
      (ac.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
  });

  await softStep('Click OK (defaults) → Scatter plot of cliffs appears', async () => {
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(45000);
    const viewers = await page.evaluate(() => Array.from(grok.shell.tv.viewers).map((v: any) => v.type));
    expect(viewers).toContain('Scatter plot');
  });

  await softStep('Click cliffs count link → cliffs grid appears', async () => {
    const hasLink = await page.evaluate(() => {
      const links = document.querySelectorAll('[name*="button-"][name*="CLIFFS"], [name*="button-"][name*="cliffs"]');
      return links.length > 0;
    });
    if (hasLink) {
      await page.evaluate(async () => {
        const btn = document.querySelector('[name*="button-"][name*="CLIFFS"], [name*="button-"][name*="cliffs"]') as HTMLElement;
        if (btn) btn.click();
        await new Promise(r => setTimeout(r, 3000));
      });
      const viewerCount = await page.evaluate(() => Array.from(grok.shell.tv.viewers).length);
      expect(viewerCount).toBeGreaterThanOrEqual(2);
    }
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
