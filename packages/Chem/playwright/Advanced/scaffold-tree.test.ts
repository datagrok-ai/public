import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

// Smoke-scoped: run.md marks the filter/structure-filter/cloned-view/invalid-structure
// sub-scenarios as SKIP (canvas-only, not scriptable). This spec verifies the two PASS
// steps: Scaffold Tree viewer launches from the Chem menu, and the magic wand generates
// a populated scaffold tree on SPGI.csv.
test('Chem: Scaffold Tree viewer smoke', async ({page}) => {
  test.setTimeout(120_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    grok.shell.addTableView(df);
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
  await waitForMolecule(page);
  await waitForChemMenu(page);

  await softStep('Scaffold Tree viewer launches from Chem menu', async () => {
    await page.evaluate(() => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.waitForFunction(() => Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .some(m => m.textContent!.trim() === 'Scaffold Tree'), null, {timeout: 10_000});
    await page.evaluate(() => {
      const st = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Scaffold Tree') as HTMLElement;
      (st.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.waitForFunction(() => Array.from((window as any).grok.shell.tv.viewers)
      .some((v: any) => /scaffold/i.test(v.type || '')), null, {timeout: 15_000});
    const hasScaffold = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some((v: any) => /scaffold/i.test(v.type || '')));
    expect(hasScaffold).toBe(true);
  });

  await softStep('Generate scaffold tree (magic wand) → nodes appear', async () => {
    await page.waitForFunction(() => {
      const st = Array.from((window as any).grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      return !!(st as any)?.root?.querySelector('.fa-magic, [title*="Generate" i]');
    }, null, {timeout: 30_000});
    const started = await page.evaluate(() => {
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      const wand = (st as any)?.root?.querySelector('.fa-magic, [title*="Generate" i]');
      if (wand) { (wand as HTMLElement).click(); return true; }
      return false;
    });
    expect(started, 'magic-wand icon present on Scaffold Tree viewer').toBe(true);
    await page.waitForFunction(() => {
      const st = Array.from((window as any).grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      return ((st as any)?.root?.querySelectorAll('.d4-tree-view-node, .d4-scaffold-tree-node').length || 0) >= 1;
    }, null, {timeout: 60_000});
    const nodeCount = await page.evaluate(() => {
      const st = Array.from(grok.shell.tv.viewers).find((v: any) => /scaffold/i.test(v.type || ''));
      return (st as any)?.root?.querySelectorAll('.d4-tree-view-node, .d4-scaffold-tree-node').length || 0;
    });
    expect(nodeCount, `scaffold tree nodes: ${nodeCount}`).toBeGreaterThanOrEqual(1);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
