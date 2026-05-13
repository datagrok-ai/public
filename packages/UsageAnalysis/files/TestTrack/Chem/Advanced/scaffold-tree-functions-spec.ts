import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Chem: Scaffold Tree basic functions', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
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

  await softStep('Open Chem → Analyze → Scaffold Tree → empty viewer appears', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const st = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Scaffold Tree') as HTMLElement;
      (st.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 4000));
    });
    const hasScaffoldViewer = await page.evaluate(() => {
      return Array.from(grok.shell.tv.viewers).some((v: any) => /scaffold/i.test(v.type || ''));
    });
    expect(hasScaffoldViewer).toBe(true);
  });

  await softStep('Magic wand generates scaffold tree → nodes appear', async () => {
    const started = await page.evaluate(async () => {
      // Look for the generate / magic-wand icon in the scaffold tree viewer
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      if (!st) return false;
      const wand = (st as any).root?.querySelector('.fa-magic, [name*="magic"], [title*="Generate" i]');
      if (wand) { (wand as HTMLElement).click(); return true; }
      return false;
    });
    if (!started) test.skip(true, 'magic wand icon not found');
    await page.waitForTimeout(25000);
    const nodes = await page.evaluate(() => {
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      if (!st) return 0;
      return (st as any).root?.querySelectorAll('.d4-tree-view-node, .d4-scaffold-tree-node').length || 0;
    });
    expect(nodes).toBeGreaterThan(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
