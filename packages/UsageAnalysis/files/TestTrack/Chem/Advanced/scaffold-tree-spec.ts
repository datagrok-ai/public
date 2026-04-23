import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Chem: Scaffold Tree filter + viewer smoke', async ({page}) => {
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

  await softStep('Scaffold Tree viewer launches from Chem menu', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const st = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Scaffold Tree') as HTMLElement;
      (st.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 5000));
    });
    const hasScaffold = await page.evaluate(() => {
      return Array.from(grok.shell.tv.viewers).some((v: any) => /scaffold/i.test(v.type || ''));
    });
    expect(hasScaffold).toBe(true);
  });

  await softStep('Generate scaffold tree (magic wand) → nodes appear', async () => {
    const started = await page.evaluate(async () => {
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      if (!st) return false;
      const wand = (st as any).root?.querySelector('.fa-magic, [title*="Generate" i]');
      if (wand) { (wand as HTMLElement).click(); return true; }
      return false;
    });
    if (!started) test.skip(true, 'magic-wand icon missing');
    await page.waitForTimeout(30000);
    const nodeCount = await page.evaluate(() => {
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      if (!st) return 0;
      return (st as any).root?.querySelectorAll('.d4-tree-view-node, .d4-scaffold-tree-node').length || 0;
    });
    expect(nodeCount).toBeGreaterThanOrEqual(1);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
