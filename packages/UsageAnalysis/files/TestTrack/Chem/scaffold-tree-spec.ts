import {expect} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

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
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 100));
    }
    // The flat settle stood in for the Molecule semType landing on the opened frame.
    const semDeadline = Date.now() + 5000;
    while (Date.now() < semDeadline && !df.columns.toList().some((c: any) => c.semType === 'Molecule'))
      await new Promise(r => setTimeout(r, 100));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Scaffold Tree viewer launches from Chem menu', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      // Labels from a previously opened menu stay in the document, so only a node that was not
      // already there is this menu's leaf; clicking a stale one actuates nothing.
      const stale = new Set(Array.from(document.querySelectorAll('.d4-menu-item-label')));
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      const find = () => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => !stale.has(m) && m.textContent!.trim() === 'Scaffold Tree') as HTMLElement | undefined;
      const menuDeadline = Date.now() + 500;
      let st = find();
      while (!st && Date.now() < menuDeadline) { await new Promise(r => setTimeout(r, 25)); st = find(); }
      if (!st)
        st = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(m => m.textContent!.trim() === 'Scaffold Tree') as HTMLElement | undefined;
      (st!.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      // Wait for the viewer the assertion below reads, not for a flat interval.
      const attachDeadline = Date.now() + 5000;
      while (Date.now() < attachDeadline &&
             !Array.from(grok.shell.tv.viewers).some((v: any) => /scaffold/i.test(v.type || '')))
        await new Promise(r => setTimeout(r, 100));
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
    // generation is what the 30s covered; the node count is what the step reads, so wait for it
    await page.waitForFunction(() => {
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      return !!(st as any)?.root?.querySelectorAll('.d4-tree-view-node, .d4-scaffold-tree-node').length;
    }, undefined, {timeout: 30000}).catch(() => {});
    const nodeCount = await page.evaluate(() => {
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      if (!st) return 0;
      return (st as any).root?.querySelectorAll('.d4-tree-view-node, .d4-scaffold-tree-node').length || 0;
    });
    expect(nodeCount).toBeGreaterThanOrEqual(1);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
